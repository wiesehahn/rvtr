# rvt_data_lgln(). Location handling and the cost guard run offline; everything
# that reads data needs the LGLN catalogue and skips itself when offline.

lgln_info <- function(p) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  list(nx = nx, ny = ny, nb = ds$getRasterCount(), gt = ds$getGeoTransform(),
       bbox = ds$bbox(), epsg = gdalraster::srs_find_epsg(ds$getProjection()),
       v = as.numeric(ds$read(1, 0, 0, nx, ny, nx, ny)),
       nodata = ds$getNoDataValue(1))
}

test_that("a point resolves to the 1 km tile it lies in, however it is given", {
  tile <- c(565000, 5720000, 566000, 5721000)
  expect_equal(rvtr:::.lgln_area(c(9.9464, 51.6317), 1)$extent, tile)
  skip_if_not_installed("sf")
  pt <- sf::st_sfc(sf::st_point(c(565400, 5720700)), crs = 25832)
  expect_equal(rvtr:::.lgln_area(pt, 1)$extent, tile)
})

test_that("an extent returns exactly itself, snapped outward to the grid", {
  skip_if_not_installed("sf")
  bb <- sf::st_bbox(c(xmin = 565200.3, ymin = 5720300, xmax = 566800.4,
                      ymax = 5721100), crs = 25832)
  expect_equal(rvtr:::.lgln_area(bb, 1)$extent, c(565200, 5720300, 566801, 5721100))
  # not widened to whole tiles
  expect_equal(rvtr:::.lgln_area(bb, 0.2)$extent, c(565200.2, 5720300, 566800.4, 5721100))
})

test_that("tile extents are read from item ids", {
  expect_equal(rvtr:::.lgln_tile_extent("dop20rgbi_32_564_5720_2_ni_2025-03-08"),
               c(564000, 5720000, 566000, 5722000))
  expect_equal(rvtr:::.lgln_tile_extent("dgm1_32_565_5720_1_ni_2016"),
               c(565000, 5720000, 566000, 5721000))
})

test_that("bad arguments and oversized requests stop before any network use", {
  expect_error(rvt_data_lgln(c(9.9, 51.6), "lidar"), "must be one of")
  expect_error(rvt_data_lgln(c(565400, 5720700)), "sf point")
  # all of Lower Saxony at 100 m: the terrain floor makes this tens of GB
  expect_error(rvt_data_lgln(c(6.6, 51.3, 11.6, 53.9), "dtm", res = 100),
               "no overview coarser than 2 m")
})

test_that("the lazy raster matches the downloaded copy", {
  skip_if_offline()
  pt <- c(9.9464, 51.6317)
  lazy <- lgln_info(rvt_data_lgln(pt, "dtm"))
  local <- lgln_info(rvt_data_lgln(pt, "dtm", download = TRUE))
  expect_equal(c(lazy$nx, lazy$ny), c(1000, 1000))
  expect_equal(lazy$epsg, "EPSG:25832")
  expect_equal(lazy$bbox, c(565000, 5720000, 566000, 5721000))
  expect_identical(lazy$v, local$v)

  # and when resampling, through the warped virtual raster
  lazy2 <- lgln_info(rvt_data_lgln(pt, "rgb", res = 1))
  local2 <- lgln_info(rvt_data_lgln(pt, "rgb", res = 1, download = TRUE))
  expect_identical(lazy2$v, local2$v)
})

test_that("a second download is served from the cache", {
  skip_if_offline()
  pt <- c(9.9464, 51.6317)
  p <- rvt_data_lgln(pt, "dtm", download = TRUE)
  before <- fs::file_info(p)$modification_time
  expect_equal(rvt_data_lgln(pt, "dtm", download = TRUE), p)
  expect_equal(fs::file_info(p)$modification_time, before)
})

test_that("the native orthophoto keeps its 0.2 m pixels on the tile", {
  skip_if_offline()
  i <- lgln_info(rvt_data_lgln(c(9.9464, 51.6317), "rgb", year = 2025))
  expect_equal(c(i$nx, i$ny, i$nb), c(5000, 5000, 3))
  expect_equal(i$bbox, c(565000, 5720000, 566000, 5721000))
})

test_that("an extent across tile borders comes back whole and seamless", {
  skip_if_offline(); skip_if_not_installed("sf")
  bb <- sf::st_bbox(c(xmin = 565500, ymin = 5720500, xmax = 566500,
                      ymax = 5721500), crs = 25832)
  r <- lgln_info(rvt_data_lgln(bb, "dtm"))
  expect_equal(c(r$nx, r$ny), c(1000, 1000))
  expect_false(any(is.na(r$v) | r$v == r$nodata))
})

test_that("lazy resampled reads are tile-invariant", {
  skip_if_offline()
  dem <- rvt_data_lgln(c(9.9464, 51.6317), "dtm", res = 2)
  whole <- lgln_info(rvt_svf(dem, reach = 20))$v
  pieces <- lgln_info(rvt_svf(dem, reach = 20, tile_size = 137))$v
  expect_identical(whole, pieces)
})

test_that("products fetched for one place at one res share a grid", {
  skip_if_offline()
  pt <- c(9.9464, 51.6317)
  dtm <- rvt_data_lgln(pt, "dtm", res = 1)
  dsm <- rvt_data_lgln(pt, "dsm", res = 1)
  rgb <- rvt_data_lgln(pt, "rgb", res = 1)
  a <- lgln_info(dtm); b <- lgln_info(dsm); c <- lgln_info(rgb)
  expect_equal(b$gt, a$gt); expect_equal(c$gt, a$gt)
  expect_equal(c(c$nx, c$ny, c$nb), c(a$nx, a$ny, 3))
  hs <- rvt_hillshade(dtm)
  on.exit(unlink(hs), add = TRUE)
  expect_equal(lgln_info(rvt_render(rvt_blend(rgb, hs, "multiply", 0.6)))$nb, 3L)
})

test_that("years are listed and selectable", {
  skip_if_offline()
  pt <- c(9.9464, 51.6317)
  y <- rvt_data_lgln_years(pt, "rgb")
  expect_true(all(c(2013, 2016, 2019, 2022, 2025) %in% y$year))
  expect_true(all(y$covers))
  old <- lgln_info(rvt_data_lgln(pt, "rgb", year = 2013, res = 10))
  new <- lgln_info(rvt_data_lgln(pt, "rgb", year = 2025, res = 10))
  expect_equal(old$gt, new$gt)
  expect_false(identical(old$v, new$v))
  expect_error(rvt_data_lgln(pt, "rgb", year = 2001), "Years that do: 2013")
})

test_that("places outside Lower Saxony are refused", {
  skip_if_offline()
  expect_error(rvt_data_lgln(c(11.58, 48.14), "dtm"), "only covers Lower Saxony")
})
