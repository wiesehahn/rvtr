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
  bdom <- rvt_data_lgln(pt, "bdom", res = 1)
  a <- lgln_info(dtm); b <- lgln_info(dsm); c <- lgln_info(rgb)
  d <- lgln_info(bdom)
  expect_equal(b$gt, a$gt); expect_equal(c$gt, a$gt); expect_equal(d$gt, a$gt)
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

test_that("a raster `x` sets the grid, and the product default doesn't fight it", {
  # A DTM rasterised from a LAS at 0.25 m: its own grid is the answer, and
  # asking for the orthophoto must not compare 0.25 against the product's
  # native 0.2. That used to make rvt_data_lgln_years(dtm, "rgb") and
  # rvt_data_lgln(dtm, "rgb") errors, though nothing had asked for a resolution.
  p <- fs::file_temp(ext = "tif")
  on.exit(unlink(p), add = TRUE)
  gdalraster::create("GTiff", p, 200, 200, 1, "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(564500, 0.25, 0, 5720550, 0, -0.25))
  ds$setProjection(gdalraster::srs_to_wkt("EPSG:25832"))
  ds$close()

  a <- rvtr:::.lgln_area(p, NULL, 0.2)          # nothing asked, rgb default
  expect_equal(a$res, 0.25)
  expect_equal(a$kind, "raster")
  expect_equal(a$extent, c(564500, 5720500, 564550, 5720550))

  # an explicit, conflicting request is still refused
  expect_error(rvtr:::.lgln_area(p, 1, 0.2), "can't differ")
  expect_equal(rvtr:::.lgln_area(p, 0.25, 0.2)$res, 0.25)

  # and the default still applies where the location brings no grid
  expect_equal(rvtr:::.lgln_area(c(9.9464, 51.6317), NULL, 0.2)$res, 0.2)
})

test_that("bdom is costed from its overview ladder, and the guard says so", {
  # bDOM is uncompressed 0.2 m Float32, so a read costs exactly the pixels GDAL
  # takes from one level times four bytes: 100 MB per km² at native, and one
  # step down per overview. Nothing else here is remotely this heavy.
  expect_equal(rvtr:::.lgln_mb_per_km2("bdom", 0.2), 100)
  expect_equal(rvtr:::.lgln_mb_per_km2("bdom", 0.4), 25)
  expect_equal(rvtr:::.lgln_mb_per_km2("bdom", 1), 6.25)
  expect_equal(rvtr:::.lgln_mb_per_km2("bdom", 2), 1.56)
  expect_equal(rvtr:::.lgln_mb_per_km2("bdom", 50), 0.39)
  # never cheaper for a finer request
  mb <- vapply(c(0.1, 0.2, 0.3, 0.5, 1, 3, 10), function(r)
    rvtr:::.lgln_mb_per_km2("bdom", r), numeric(1))
  expect_false(is.unsorted(rev(mb)))

  # the guard fires before any network use, naming `res` rather than the
  # terrain hint, which belongs to dtm/dsm
  expect_error(rvt_data_lgln(c(9.90, 51.60, 9.96, 51.64), "bdom"),
               "20 times the orthophoto")
  expect_error(rvt_data_lgln(c(6.6, 51.3, 11.6, 53.9), "dtm", res = 100),
               "no overview coarser than 2 m")
})

test_that("bdom item ids parse as 1 km tiles and the product is accepted", {
  expect_equal(rvtr:::.lgln_tile_extent("bdom20_32_565_5720_1_ni_2025-03-08"),
               c(565000, 5720000, 566000, 5721000))
  expect_error(rvt_data_lgln(c(9.9, 51.6), "lidar"), '"bdom"')
  expect_error(rvt_data_lgln_years(c(9.9, 51.6), "lidar"), '"bdom"')
})

test_that("the image-based surface model comes back on the shared grid", {
  skip_if_offline()
  pt <- c(9.9464, 51.6317)
  # res = 1 on purpose: native 0.2 m would pull 100 MB for this one tile
  b <- lgln_info(rvt_data_lgln(pt, "bdom", res = 1))
  d <- lgln_info(rvt_data_lgln(pt, "dtm"))
  expect_equal(c(b$nx, b$ny), c(d$nx, d$ny))
  expect_equal(b$gt, d$gt)
  expect_equal(b$epsg, "EPSG:25832")
  # a surface model, so it stands above the terrain somewhere
  expect_gt(max(b$v - d$v, na.rm = TRUE), 5)

  y <- rvt_data_lgln_years(pt, "bdom")
  expect_true(all(c(2022, 2025) %in% y$year))
  expect_true(all(y$covers))
  expect_true(all(y$year >= 2021))            # bDOM starts in 2021
  # how well it describes a canopy depends on the flight season, so nothing
  # here pins that - see ?rvt_data_lgln and the measurements in CLAUDE.md
})
