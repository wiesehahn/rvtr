# terra SpatRasters as input. Everything funnels through .as_path(), so these
# check the conversion rule rather than every function.

dem <- system.file("extdata", "dtm1.tif", package = "rvtr")

read_all <- function(p, band = 1L) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  as.numeric(ds$read(band, 0, 0, nx, ny, nx, ny))
}

test_that("a file-backed SpatRaster hands over its own file, with no copy", {
  skip_if_not_installed("terra")
  r <- terra::rast(dem)
  expect_equal(as.character(rvtr:::.as_path(r)), as.character(rvtr:::.as_path(dem)))

  # and the result is identical to computing from the path
  expect_identical(read_all(rvt_slope(r)), read_all(rvt_slope(dem)))
})

test_that("a SpatRaster computed in memory is written out and used", {
  skip_if_not_installed("terra")
  r <- terra::rast(dem) + 100                       # in memory, no source file
  expect_true(terra::inMemory(r))

  p <- rvtr:::.as_path(r)
  on.exit(unlink(p), add = TRUE)
  expect_false(as.character(p) == as.character(dem))
  expect_equal(read_all(p), read_all(dem) + 100, tolerance = 1e-6)

  # slope is unchanged by adding a constant, so it must match the original.
  # Not exactly: terra holds the computed raster in double and writes Float64,
  # where the source file is Float32.
  expect_equal(read_all(rvt_slope(r)), read_all(rvt_slope(dem)), tolerance = 1e-3)
})

test_that("a band subset keeps its own band, not band 1 of the file", {
  skip_if_not_installed("terra")
  # a 3-band file whose band 2 is the DEM raised by 100
  multi <- tempfile(fileext = ".tif")
  terra::writeRaster(c(terra::rast(dem), terra::rast(dem) + 100,
                       terra::rast(dem) + 200), multi)
  on.exit(unlink(multi), add = TRUE)

  sub <- terra::rast(multi)[[2]]
  p <- rvtr:::.as_path(sub)
  on.exit(unlink(p), add = TRUE)
  expect_false(as.character(p) == as.character(multi))     # not the file itself
  expect_equal(read_all(p), read_all(multi, band = 2), tolerance = 1e-6)

  # the whole file, by contrast, is passed straight through
  expect_equal(as.character(rvtr:::.as_path(terra::rast(multi))),
               as.character(rvtr:::.as_path(multi)))
})

test_that("a cropped SpatRaster keeps its own extent", {
  skip_if_not_installed("terra")
  r <- terra::rast(dem)
  e <- as.vector(terra::ext(r))
  crp <- terra::crop(r, terra::ext(e[1] + 100, e[2] - 100, e[3] + 100, e[4] - 100))
  p <- rvtr:::.as_path(crp)
  on.exit(unlink(p), add = TRUE)

  i <- rvtr:::.dem_info(p)
  expect_equal(c(i$nx, i$ny), c(800, 800))
  expect_equal(unname(as.vector(terra::ext(terra::rast(p)))),
               unname(e + c(100, -100, 100, -100)))
})

test_that("SpatRasters work through the blending and plotting entry points", {
  skip_if_not_installed("terra")
  r <- terra::rast(dem)
  hs <- rvt_hillshade(r)
  svf <- rvt_svf(r, reach = 10)
  out <- rvt_render(rvt_blend(terra::rast(hs), terra::rast(svf), "multiply", opacity = 0.5))
  expect_identical(read_all(out),
                   read_all(rvt_render(rvt_blend(hs, svf, "multiply", opacity = 0.5))))
})

test_that("a SpatRaster of any data type can be materialised", {
  skip_if_not_installed("terra")
  # A masked orthophoto - terra writes it as an integer type, and the scratch
  # copy used to ask GDAL for PREDICTOR=3, which is Float32/64 only:
  # "[writeRaster] failed writing GTiff file". terra::datatype() is "" for an
  # in-memory raster, so the type cannot be checked before writing; the scratch
  # is uncompressed instead.
  d <- terra::rast(dem)
  rgb <- terra::rast(terra::ext(d), nrows = 50, ncols = 50, nlyrs = 3,
                     crs = terra::crs(d), vals = rep(0:255, length.out = 7500))
  masked <- terra::mask(rgb, rgb[[1]] > 10, maskvalue = FALSE)

  p <- rvtr:::.as_path(masked)
  on.exit(unlink(p), add = TRUE)
  expect_true(file.exists(p))
  expect_equal(rvtr:::.nbands(p), 3L)

  img <- rvt_image(masked, fs::file_temp(ext = "webp"), range = c(0, 255))
  expect_true(file.exists(img))
})
