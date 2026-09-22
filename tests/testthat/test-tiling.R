# Structural invariants that hold regardless of any reference implementation:
# a mosaic of several files equals the merged raster bit for bit, tiling
# changes nothing, resampling only coarsens, and every distance parameter is
# in map units rather than pixels.
#
# These are the guarantees the package is built on, so they are checked
# separately for the horizon-search kernel and the integral-image kernel -
# different code paths with different overlap logic.

test_that("a mosaic gives bit-identical results to the whole raster", {
  # the tile cut into quadrants (quad_rasters(), helper-rvtr.R), then
  # processed as a mosaic
  expect_identical(read_band(rvt_svf(quad_rasters())), read_band(fx_svf()))
})

test_that("tiling does not change the result", {
  # a tile size that forces several uneven tiles plus overlap handling
  expect_equal(read_band(rvt_svf(dem, tile_size = 137)),
               read_band(fx_svf()),
               tolerance = 1e-12)
})

test_that("tiling does not change the MSRM result", {
  # separate from the horizon-kernel check above: MSRM uses a different
  # kernel (integral image, not horizon search), worth its own invariant
  expect_equal(read_band(rvt_msrm(dem, tile_size = 137)),
               read_band(fx_msrm()),
               tolerance = 1e-12)
})

test_that("rvt_resample only coarsens, and coarsens correctly", {
  p <- rvt_resample(dem, 4)
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  nx <- ds$getRasterXSize(); r <- ds$res()
  ds$close()
  expect_equal(nx, 250L)
  expect_equal(unname(r[1]), 4)

  # detail cannot be invented
  expect_error(rvt_resample(dem, 0.5), "finer than the source")
  # asking for the source resolution is allowed - a no-op copy, same grid
  same <- rvt_resample(dem, 1)
  ds <- methods::new(gdalraster::GDALRaster, same, read_only = TRUE)
  expect_equal(ds$getRasterXSize(), 1000L)
  ds$close()

  # the elevations must stay elevations - no NoData bled into the mean
  src <- read_band(dem)
  got <- read_band(p)
  expect_gte(min(got, na.rm = TRUE), min(src, na.rm = TRUE))
  expect_lte(max(got, na.rm = TRUE), max(src, na.rm = TRUE))
})

test_that("a metric on a resampled DEM is still tile-invariant", {
  # the reason rvt_resample() writes a real file instead of a lazy VRT:
  # GDAL's on-the-fly resampling can be window-dependent, which would make
  # results depend on tile_size. Both an integer and a non-integer factor.
  for (res in c(4, 2.5)) {
    p <- rvt_resample(dem, res)
    expect_identical(read_band(rvt_svf(p, tile_size = 37)),
                     read_band(rvt_svf(p)))
  }
})

test_that("categorical resampling keeps whole class codes", {
  # pins `-ovr NONE`: without it GDAL serves the source COG's AVERAGE-built
  # overviews and silently ignores method = "mode", returning class codes
  # like 4.31 that correspond to no geomorphon at all.
  g <- rvt_geomorphons(dem)
  classes <- sort(unique(as.vector(read_band(g))))
  classes <- classes[!is.na(classes)]

  m <- read_band(rvt_resample(g, 4, method = "mode"))
  got <- sort(unique(as.vector(m)))
  got <- got[!is.na(got)]
  expect_true(all(got == round(got)))
  expect_true(all(got %in% classes))
})

test_that("tiling does not change time in daylight", {
  args <- list(dem, dates = as.Date("2025-06-21"), reach = 20, time_step = 30)
  expect_equal(read_band(do.call(rvt_daylight, c(args, list(tile_size = 137)))),
               read_band(do.call(rvt_daylight, args)),
               tolerance = 1e-12)
})

test_that("tiling and mosaicking do not change a multi-resolution result", {
  one <- read_band(rvt_svf(dem, reach = 100, pyramid_px = 25))
  expect_identical(read_band(rvt_svf(dem, reach = 100, pyramid_px = 25,
                                      tile_size = 137)), one)

  # the same four quadrant files the plain mosaic test above uses, so the
  # coarse levels are built from a VRT over several files rather than one
  expect_identical(read_band(rvt_svf(quad_rasters(), reach = 100,
                                      pyramid_px = 25)), one)
})

test_that("distances snap to whole cells by a stated rule", {
  # halfway cases must not inherit round()'s round-half-to-even, which would
  # send 2.5 down and 3.5 up
  expect_equal(.cells(5, 2, "d", "outer"), 3L)     # 2.5 -> 3, not 2
  expect_equal(.cells(7, 2, "d", "outer"), 4L)     # 3.5 -> 4
  expect_equal(.cells(4, 2, "d", "outer"), 2L)
  expect_equal(.cells(10, 1, "d", "outer"), 10L)   # exact, and silent
  expect_silent(.cells(10, 1, "d", "outer"))

  # below one cell: an outer radius cannot be honoured at all
  expect_error(.cells(1, 2, "radius", "outer"), "smaller than one cell")
  # but an inner radius or skip may legitimately be zero
  expect_equal(.cells(1, 2, "skip", "inner"), 1L)
  expect_equal(.cells(0.2, 2, "skip", "inner"), 0L)
  # and a step finer than the grid just means every cell
  expect_equal(.cells(0.2, 2, "step", "step"), 1L)

  # say so when the grid could not honour the request closely
  expect_message(.cells(5, 2, "radius", "outer"), "snapped to 3 cells")
})

test_that("metric parameters are distances, not pixel counts", {
  # a 2 m raster: the same request must mean the same ground distance
  n <- 120L
  m <- matrix(rnorm(n * n, 100, 2), n, n)
  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = n,
                           ysize = n, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(0, 2, 0, 0, 0, -2))        # 2 m cells
  ds$write(1L, 0, 0, n, n, as.vector(t(m)))
  ds$close()
  on.exit(unlink(p), add = TRUE)

  # asking for less than one cell is refused rather than silently rounded
  expect_error(rvt_svf(p, reach = 1), "smaller than one cell")
  expect_error(rvt_slrm(p, radius = 1), "smaller than one cell")
  # and a request the grid cannot hit exactly says so
  expect_message(rvt_svf(p, reach = 5), "snapped")

  # 20 m of smoothing on a 2 m raster must equal 10 cells, which is what the
  # same 20 m means on the 1 m sample tile
  expect_equal(read_band(rvt_slrm(p, radius = 20)),
               read_band(rvt_slrm(p, radius = 20.4)), tolerance = 1e-12)
})

test_that("tiling does not change a multi-resolution daylight or shadow", {
  one <- read_band(rvt_shadow(dem, reach = 100, pyramid_px = 25))
  expect_identical(read_band(rvt_shadow(dem, reach = 100, pyramid_px = 25,
                                         tile_size = 137)), one)

  args <- list(dem, dates = as.Date("2025-06-21"), time_step = 60,
               num_directions = 8, reach = 100, pyramid_px = 25)
  expect_identical(read_band(do.call(rvt_daylight, c(args, list(tile_size = 137)))),
                   read_band(do.call(rvt_daylight, args)))

  si <- list(dem, num_directions = 8, reach = 100, pyramid_px = 25)
  expect_identical(read_band(do.call(rvt_sky_illumination, c(si, list(tile_size = 137)))),
                   read_band(do.call(rvt_sky_illumination, si)))
})

test_that("max resampling really takes the block maximum", {
  # gdal_translate cannot do -r max and silently substitutes nearest, which
  # would make the pyramid's coarse levels miss narrow distant obstructions
  # entirely; rvt_resample routes these methods through gdalwarp instead.
  z <- matrix(0, 64, 64)
  z[30:31, ] <- 50                      # a 2-cell ridge, thinner than one block
  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = 64,
                           ysize = 64, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(557000, 1, 0, 5700000, 0, -1))
  ds$setProjection(gdalraster::srs_to_wkt("EPSG:25832"))
  ds$write(1L, 0, 0, 64, 64, as.vector(t(z)))
  ds$close()
  on.exit(unlink(p), add = TRUE)

  got <- read_band(rvt_resample(p, 8, method = "max"))
  expect_equal(max(got), 50)                  # the ridge survives
  expect_equal(sum(got == 50), 8)             # in exactly one row of blocks
  # nearest, by contrast, steps straight over it
  expect_equal(max(read_band(rvt_resample(p, 8, method = "near"))), 0)

  # and min is the mirror image, which is what negative openness needs
  expect_equal(max(read_band(rvt_resample(p, 8, method = "min"))), 0)
})

