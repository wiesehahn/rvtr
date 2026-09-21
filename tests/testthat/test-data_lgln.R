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

test_that("coverage is the exact area a year's tiles fill", {
  # tiles as .lgln_search() returns them, on the whole-kilometre grid
  tiles <- function(...) {
    e <- do.call(rbind, list(...))
    data.frame(xmin = e[, 1], ymin = e[, 2], xmax = e[, 3], ymax = e[, 4])
  }
  one <- tiles(c(565000, 5720000, 566000, 5721000))
  ext <- c(565000, 5720000, 566000, 5721000)

  expect_equal(rvtr:::.lgln_coverage(one, ext), 1)
  expect_true(rvtr:::.lgln_covers(one, ext))

  # an extent reaching one tile east: half of it has data
  expect_equal(rvtr:::.lgln_coverage(one, c(565000, 5720000, 567000, 5721000)), 0.5)
  expect_false(rvtr:::.lgln_covers(one, c(565000, 5720000, 567000, 5721000)))

  # a quarter, reaching east and north
  expect_equal(rvtr:::.lgln_coverage(one, c(565000, 5720000, 567000, 5722000)), 0.25)

  # two tiles side by side fill the pair
  two <- tiles(c(565000, 5720000, 566000, 5721000),
               c(566000, 5720000, 567000, 5721000))
  expect_equal(rvtr:::.lgln_coverage(two, c(565000, 5720000, 567000, 5721000)), 1)

  # an L leaves the fourth quadrant empty
  ell <- tiles(c(565000, 5720000, 566000, 5721000),
               c(566000, 5720000, 567000, 5721000),
               c(565000, 5721000, 566000, 5722000))
  expect_equal(rvtr:::.lgln_coverage(ell, c(565000, 5720000, 567000, 5722000)), 0.75)

  # overlapping tiles are counted once, not summed - the probe-free sweep is
  # what makes that free
  ovl <- tiles(c(565000, 5720000, 566000, 5721000),
               c(565500, 5720000, 566500, 5721000))
  expect_equal(rvtr:::.lgln_coverage(ovl, c(565000, 5720000, 566500, 5721000)), 1)

  # an extent well inside one tile, and one entirely outside
  expect_equal(rvtr:::.lgln_coverage(one, c(565100, 5720100, 565200, 5720200)), 1)
  expect_equal(rvtr:::.lgln_coverage(one, c(570000, 5730000, 571000, 5731000)), 0)
  expect_equal(rvtr:::.lgln_coverage(NULL, ext), 0)

  # a gap narrower than the old 500 m probe step, which used to read as full
  gap <- tiles(c(565000, 5720000, 565400, 5721000),
               c(565500, 5720000, 566000, 5721000))
  expect_equal(rvtr:::.lgln_coverage(gap, ext), 0.9)

  # the same tile under two flight dates is one tile
  dup <- tiles(c(565000, 5720000, 566000, 5721000),
               c(565000, 5720000, 566000, 5721000))
  expect_equal(rvtr:::.lgln_coverage(dup, c(565000, 5720000, 567000, 5721000)), 0.5)
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
  expect_error(rvt_data_lgln(pt, "rgb", year = 2001), "Years present: 2013")
  expect_true(all(y$coverage == 1))
})

test_that("a named year is used on its own, hole and all", {
  skip_if_offline()
  skip_if_not_installed("sf")
  # A real seam between flight campaigns, west of Celle: *no* single year
  # covers this 2 x 1 km box - three surface-model years each reach half of
  # it - which is exactly the case the old code refused outright.
  bb <- sf::st_bbox(c(xmin = 583000, ymin = 5808000, xmax = 585000,
                      ymax = 5809000), crs = 25832)
  y <- rvt_data_lgln_years(bb, "dsm")
  skip_if(all(y$covers), "this location is no longer a seam")
  yr <- y$year[y$coverage > 0.05 & y$coverage < 0.95][1]

  expect_message(got <- rvt_data_lgln(bb, "dsm", year = yr, res = 10), "covers")
  v <- lgln_info(got)$v
  expect_true(anyNA(v) && any(!is.na(v)))
  # the hole is exactly the share the catalogue promised
  expect_equal(mean(!is.na(v)), y$coverage[y$year == yr], tolerance = 0.02)

  # and partial = FALSE refuses it, naming both ways out
  expect_error(rvt_data_lgln(bb, "dsm", year = yr, res = 10, partial = FALSE),
               "partial = TRUE")
})

test_that("year = NULL gives each tile its latest flight, across years", {
  skip_if_offline()
  skip_if_not_installed("sf")
  # The same seam. No year covers it alone, but together they do, so the
  # default must return a complete raster assembled from several years -
  # newest where each tile was last flown.
  bb <- sf::st_bbox(c(xmin = 583000, ymin = 5808000, xmax = 585000,
                      ymax = 5809000), crs = 25832)
  y <- rvt_data_lgln_years(bb, "dsm")
  skip_if(all(y$covers), "this location is no longer a seam")

  expect_message(got <- rvt_data_lgln(bb, "dsm", res = 10), "combines 2 years")
  expect_false(anyNA(lgln_info(got)$v))
  # complete, so partial = FALSE must accept it even though no single year does
  expect_no_error(rvt_data_lgln(bb, "dsm", res = 10, partial = FALSE))

  # each tile takes the newest year that holds it, not the newest year overall
  items <- rvtr:::.lgln_search("dsm", as.numeric(sf::st_bbox(bb)))
  picked <- rvtr:::.lgln_pick(items, as.numeric(sf::st_bbox(bb)), NULL, "dsm")
  expect_length(unique(picked$items$year), 2L)
  expect_equal(picked$coverage, 1)
  expect_equal(picked$year, "2016-2019")
  for (i in seq_len(nrow(picked$items))) {
    same <- items[items$xmin == picked$items$xmin[i] &
                    items$ymin == picked$items$ymin[i], , drop = FALSE]
    expect_equal(picked$items$date[i], max(same$date))
  }
})

test_that("where one year covers everything, that year is what you get", {
  skip_if_offline()
  pt <- c(9.9464, 51.6317)
  y <- rvt_data_lgln_years(pt, "rgb")
  a <- lgln_info(rvt_data_lgln(pt, "rgb", res = 10))
  b <- lgln_info(rvt_data_lgln(pt, "rgb", year = max(y$year), res = 10))
  expect_equal(a$v, b$v)
  # fully covered by one year, so neither the coverage nor the mixed-year
  # message fires (the licence one can, hence matching text not expect_silent)
  msgs <- testthat::capture_messages(
    rvt_data_lgln(pt, "rgb", res = 10, partial = FALSE))
  expect_false(any(grepl("covers|combines", msgs)))
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
