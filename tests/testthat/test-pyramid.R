# The multi-resolution horizon search: near terrain at full resolution,
# distant terrain from coarser copies of the same DEM (src/levels.h,
# R/pyramid.R). Compared against an exact full-resolution search of the same
# reach, which is what a pyramid_px covering the whole reach in one level
# gives.
#
# Always away from the raster edge: within `reach` of the boundary both
# methods are extrapolating terrain that is not there, and they do it
# differently (level 0 reflects, coarse levels replicate).

test_that("the pyramid plan covers the requested reach", {
  p <- .pyramid_plan(xres = 1, reach = 400, px = 100, factor = 4)
  expect_equal(length(p), 2L)
  expect_equal(p[[1]]$fac, 1L)          # level 0 is always the DEM itself
  expect_equal(p[[1]]$outer, 100)
  expect_equal(p[[2]]$fac, 4L)
  expect_equal(p[[2]]$outer, 400)       # reaches exactly as far as asked
  expect_equal(p[[2]]$inner, 100)       # and starts where level 0 stopped

  # levels stay contiguous and the last one always lands on `reach`
  for (cfg in list(c(0.25, 1600, 100, 4), c(1, 250, 50, 2), c(2, 5000, 100, 8))) {
    p <- .pyramid_plan(cfg[1], cfg[2], cfg[3], cfg[4])
    expect_equal(p[[length(p)]]$outer, cfg[2])
    for (i in seq_along(p)[-1]) expect_equal(p[[i]]$inner, p[[i - 1]]$outer)
  }
  expect_error(.pyramid_plan(1, -5, 100, 4), "positive")
  expect_error(.pyramid_plan(1, 400, 100, 1), "at least 2")
})

test_that("a multi-resolution search matches a full-resolution one", {
  # Compared away from the raster edge: within `reach` of the boundary both
  # methods are extrapolating terrain that isn't there, and they do it
  # differently (level 0 reflects, coarse levels replicate).
  d <- crop_raster(dem)
  interior <- function(m, k) m[(k + 1):(nrow(m) - k), (k + 1):(ncol(m) - k)]
  err <- function(a, b, k = 300) {
    a <- interior(read_band(a), k); b <- interior(read_band(b), k)
    ok <- !is.na(a) & !is.na(b)
    mean(abs(a[ok] - b[ok]))
  }
  # reach 200 with the default pyramid_px gives two levels; pyramid_px = 200
  # covers the same reach in one, which is the exact search to compare against
  expect_lt(err(rvt_svf(d, reach = 200, pyramid_px = 200),
                rvt_svf(d, reach = 200)), 1e-4)
  expect_lt(err(fx_opns_exact(), rvt_openness(d, reach = 200)), 0.02)
  # Negative openness inverts every level, which only works if the coarse
  # levels take the mirror statistic (min for max - see .mirror_method()).
  # Its tolerance is deliberately looser than positive openness's 0.02: on
  # this terrain the error is 0.020 against 0.004, and it does not shrink with
  # a wider margin, so it is not an edge effect. Inverting turns narrow
  # incised gullies into narrow ridges, and narrow features are exactly what
  # coarsening represents least well - the same asymmetry the distant-wall
  # measurements show.
  expect_lt(err(rvt_openness_negative(d, reach = 200, pyramid_px = 200),
                rvt_openness_negative(d, reach = 200)), 0.03)
})

test_that("a coarser pyramid is worse but still close", {
  # Reference and candidates must share a reach - the comparison is about
  # pyramid_px alone. 200 m keeps a far field worth measuring while costing a
  # quarter of what 400 m did.
  d <- crop_raster(dem)
  interior <- function(m, k = 300) m[(k + 1):(nrow(m) - k), (k + 1):(ncol(m) - k)]
  ref <- interior(read_band(fx_opns_exact()))
  e <- function(px) {
    got <- interior(read_band(rvt_openness(d, reach = 200, pyramid_px = px)))
    ok <- !is.na(got) & !is.na(ref)
    mean(abs(got[ok] - ref[ok]))
  }
  # pyramid_px is the accuracy knob: more full-resolution reach, less error
  expect_lt(e(100), e(50))
  expect_lt(e(50), e(25))
})

test_that("daylight, shadow and sky illumination also search multi-resolution", {
  # Each of these has its own kernel, all three now driven by the same
  # RvtLevels accumulation (src/levels.h). Compared against an exact
  # full-resolution search of the same reach - which is what a pyramid_px
  # large enough to cover the whole reach in one level gives - away from the
  # raster edge, where the two extrapolate differently.
  d <- crop_raster(dem)
  interior <- function(m, k = 300) m[(k + 1):(nrow(m) - k), (k + 1):(ncol(m) - k)]
  err <- function(a, b) {
    a <- interior(read_band(a)); b <- interior(read_band(b))
    ok <- !is.na(a) & !is.na(b)
    mean(abs(a[ok] - b[ok]))
  }

  expect_lt(err(rvt_shadow(d, reach = 200, pyramid_px = 200),
                rvt_shadow(d, reach = 200)), 0.01)

  expect_lt(err(rvt_sky_illumination(d, num_directions = 8, reach = 200,
                                      pyramid_px = 200),
                rvt_sky_illumination(d, num_directions = 8, reach = 200)), 1e-3)

  expect_lt(err(rvt_daylight(d, dates = as.Date("2025-06-21"), time_step = 60,
                              num_directions = 8, reach = 200, pyramid_px = 200),
                rvt_daylight(d, dates = as.Date("2025-06-21"), time_step = 60,
                              num_directions = 8, reach = 200)), 0.05)
})

test_that("pyramid_method chooses how the coarse levels are built", {
  # max keeps a narrow obstruction, q3 discards it - that is the whole trade,
  # and which is right depends on whether the surface is bare earth or canopy
  z <- matrix(0, 64, 64)
  z[30:31, ] <- 50
  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = 64,
                           ysize = 64, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(557000, 1, 0, 5700000, 0, -1))
  ds$setProjection(gdalraster::srs_to_wkt("EPSG:25832"))
  ds$write(1L, 0, 0, 64, 64, as.vector(t(z)))
  ds$close()
  on.exit(unlink(p), add = TRUE)
  expect_equal(max(read_band(rvt_resample(p, 8, method = "max"))), 50)
  expect_equal(max(read_band(rvt_resample(p, 8, method = "q3"))), 0)

  # the argument reaches the levels: a rough surface coarsened by q3 must see
  # less obstruction than the same surface coarsened by max. Needs a CRS,
  # since these statistics are done with gdalwarp.
  set.seed(3)
  rough <- matrix(pmax(0, rnorm(400 * 400, 0, 6)), 400, 400)
  r <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = r, xsize = 400,
                           ysize = 400, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, r, read_only = FALSE)
  ds$setGeoTransform(c(557000, 1, 0, 5700000, 0, -1))
  ds$setProjection(gdalraster::srs_to_wkt("EPSG:25832"))
  ds$write(1L, 0, 0, 400, 400, as.vector(t(rough)))
  ds$close()
  on.exit(unlink(r), add = TRUE)

  mx <- read_band(rvt_svf(r, reach = 200, pyramid_method = "max"))
  q3 <- read_band(rvt_svf(r, reach = 200, pyramid_method = "q3"))
  expect_gt(mean(q3), mean(mx))

  # negative openness inverts the terrain, so its levels take the mirror
  # statistic; a method with no mirror is refused rather than silently wrong
  expect_equal(.mirror_method("max"), "min")
  expect_equal(.mirror_method("q3"), "q1")
  expect_equal(.mirror_method("med"), "med")
  expect_equal(.mirror_method("cubic"), "cubic")
  expect_error(.mirror_method("rms"), "does not commute")
  expect_error(rvt_openness_negative(dem, reach = 400, pyramid_method = "rms"),
               "does not commute")
})

