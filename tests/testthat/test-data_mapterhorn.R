# rvt_data_mapterhorn(). Zoom arithmetic, CRS choice and argument checks run
# offline; everything that reads coverage or tiles skips itself when offline.

mt_info <- function(p) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  list(nx = nx, ny = ny, gt = ds$getGeoTransform(), bbox = ds$bbox(),
       epsg = gdalraster::srs_find_epsg(ds$getProjection()),
       v = as.numeric(ds$read(1, 0, 0, nx, ny, nx, ny)))
}

test_that("the zoom rule reproduces where Mapterhorn's tiles stop", {
  # each pair measured against the live server: tiles exist at this zoom and
  # 404 one zoom further
  expect_equal(rvtr:::.mt_zoom_for(1, 51.6317), 16L)     # Lower Saxony, 1 m
  expect_equal(rvtr:::.mt_zoom_for(1, 59.91), 16L)       # Norway, 1 m
  expect_equal(rvtr:::.mt_zoom_for(0.5, 46.02), 17L)     # Switzerland, 0.5 m
  expect_equal(rvtr:::.mt_zoom_for(0.25, 48.78), 18L)    # Baden-Wuerttemberg, 0.25 m
  expect_equal(rvtr:::.mt_zoom_for(2, -42.9), 15L)       # Tasmania, 2 m
  expect_equal(rvtr:::.mt_zoom_for(5, 52.1), 14L)        # Netherlands, 5 m
  # vectorised, and never beyond what a raster dimension can hold
  expect_equal(rvtr:::.mt_zoom_for(c(1, 30), 0), c(17L, 12L))
  expect_equal(rvtr:::.mt_zoom_for(1e-6, 0), rvtr:::.mt_max_zoom)
})

test_that("levels stack the global base under each source, capped by res", {
  src <- data.frame(source = c("fine", "glo30"), resolution = c(1, 30))
  ll <- c(9.93, 51.62, 9.96, 51.64)
  expect_equal(rvtr:::.mt_levels(src, 1, ll), c(12L, 16L))
  # a coarse request never reads finer tiles than it needs
  expect_equal(rvtr:::.mt_levels(src, 5, ll), c(12L, 14L))
  expect_equal(rvtr:::.mt_levels(src, 100, ll), 9L)
  # no sources at all (open ocean): the base alone
  expect_equal(rvtr:::.mt_levels(src[0, ], 30, ll), 11L)
})

test_that("the output CRS is the local UTM zone unless one is given", {
  expect_equal(rvtr:::.mt_utm(9.9, 51.6), "EPSG:32632")
  expect_equal(rvtr:::.mt_utm(-70.6, -33.4), "EPSG:32719")
  expect_equal(rvtr:::.mt_utm(180, 10), "EPSG:32660")
  expect_equal(rvtr:::.mt_utm(-180, 10), "EPSG:32601")
  a <- rvtr:::.mt_area(c(9.9464, 51.6317), 1000, NULL)
  expect_equal(a$crs, "EPSG:32632")
  expect_equal(a$extent[3] - a$extent[1], 1000)
  expect_true(a$point)
  expect_equal(rvtr:::.mt_area(c(9.9464, 51.6317), 1000, 25832)$crs, "EPSG:25832")
})

test_that("bad arguments stop before any network use", {
  expect_error(rvt_data_mapterhorn(c(9.9, 51.6), crs = 4326), "projected")
  expect_error(rvt_data_mapterhorn(c(9.9, 51.6), crs = "EPSG:999999"), "recognises")
  expect_error(rvt_data_mapterhorn(c(565400, 5720700)), "sf point")
  expect_error(rvt_data_mapterhorn(c(9.9, 51.6), res = -1), "`res`")
  expect_error(rvt_data_mapterhorn(c(9.9, 51.6), size = 0), "`size`")
  expect_error(rvt_data_mapterhorn(c(0, 87)), "85.05")
})

test_that("a point returns an exact square on the grid, in local UTM", {
  skip_if_offline()
  i <- mt_info(rvt_data_mapterhorn(c(9.9464, 51.6317)))
  expect_equal(c(i$nx, i$ny), c(1000, 1000))
  expect_equal(i$epsg, "EPSG:32632")
  expect_equal(i$bbox %% 1, c(0, 0, 0, 0))
  expect_false(anyNA(i$v))
})

test_that("heights agree with the survey's own 1 m terrain model", {
  skip_if_offline()
  lgln <- rvt_data_lgln(c(9.9464, 51.6317), "dtm")
  a <- mt_info(lgln)
  b <- mt_info(rvt_data_mapterhorn(lgln))       # on exactly the LGLN grid
  expect_equal(b$gt, a$gt)
  expect_equal(c(b$nx, b$ny), c(a$nx, a$ny))
  # measured: median |diff| 0.017 m, cor 0.999997
  expect_lt(median(abs(b$v - a$v)), 0.05)
  expect_gt(cor(a$v, b$v), 0.9999)

  # resampled: 0.033 m now, 0.098 m when the encoded RGB was being resampled
  a10 <- mt_info(rvt_data_lgln(c(9.9464, 51.6317), "dtm", res = 10))
  b10 <- mt_info(rvt_data_mapterhorn(rvt_data_lgln(c(9.9464, 51.6317), "dtm", res = 10)))
  expect_lt(median(abs(b10$v - a10$v)), 0.06)
})

test_that("resampling never touches the encoded RGB (no contour-line spikes)", {
  skip_if_offline()
  # Interpolating Terrarium RGB across a channel wrap (G every 1 m, R every
  # 256 m) makes spikes that draw contour lines in a hillshade, worst on steep
  # ground. The output may not have more Laplacian spikes than the decoded
  # tiles have at their own resolution: measured 3,434 against 4,301 on 6 km,
  # and 40,053 when the mosaic resampled the RGB.
  spikes <- function(i, thr = 30) {
    m <- matrix(i$v, i$ny, i$nx, byrow = TRUE)
    r <- 2:(i$ny - 1); k <- 2:(i$nx - 1)
    sum(abs(m[r - 1, k] + m[r + 1, k] + m[r, k - 1] + m[r, k + 1] - 4 * m[r, k]) > thr)
  }
  pt <- c(7.6586, 45.9763)                        # the Matterhorn
  p <- rvtr:::.mt_plan(pt, 10, 3000, NULL)
  z <- max(rvtr:::.mt_levels(p$sources, p$res, p$ll))
  merc <- gdalraster::bbox_transform(p$extent, p$crs, "EPSG:3857")
  native <- tempfile(fileext = ".tif")
  on.exit(unlink(native), add = TRUE)
  gdalraster::translate(rvtr:::.mt_level_vrt(z), native, quiet = TRUE,
                        cl_arg = c("-projwin", merc[1], merc[4], merc[3], merc[2]))
  out <- mt_info(rvt_data_mapterhorn(pt, res = 10, size = 3000))
  expect_lte(spikes(out), spikes(mt_info(native)))
})

test_that("where finer data stops, coarser data fills in with no holes", {
  skip_if_offline()
  # from Germany (1 m) 13 km into the Netherlands, where only the 30 m base
  # exists: the finest zoom alone leaves 45% of this NoData
  bb <- c(6.83, 52.245, 7.08, 52.26)
  p <- rvtr:::.mt_plan(bb, 10, 1000, NULL)
  lv <- rvtr:::.mt_levels(p$sources, p$res, p$ll)
  expect_gt(length(lv), 1L)
  expect_true(anyNA(mt_info(rvtr:::.mt_vrt(max(lv), p$extent, p$crs, p$res))$v))
  expect_false(anyNA(mt_info(rvt_data_mapterhorn(bb, res = 10))$v))
})

test_that("open ocean has no sources but still returns the base", {
  skip_if_offline()
  expect_equal(nrow(rvt_data_mapterhorn_sources(c(-30, 40))), 0L)
  i <- mt_info(rvt_data_mapterhorn(c(-30, 40), res = 100, size = 5000))
  expect_equal(c(i$nx, i$ny), c(50, 50))
  expect_false(anyNA(i$v))
})

test_that("sources are listed finest first, and res cannot be finer", {
  skip_if_offline()
  s <- rvt_data_mapterhorn_sources(c(9.9464, 51.6317))
  expect_equal(s$source, c("deniedersachsen", "glo30"))
  expect_equal(s$resolution, c(1, 30))
  expect_error(rvt_data_mapterhorn(c(9.9464, 51.6317), res = 0.5), "only\\s+interpolate")
})

test_that("oversized requests stop before reading any tiles", {
  skip_if_offline()
  expect_error(rvt_data_mapterhorn(c(9.5, 51.5, 10.5, 52.0)), "max_mb")
})

test_that("the lazy raster matches the downloaded copy, and reads tile-invariantly", {
  skip_if_offline()
  pt <- c(9.9464, 51.6317)
  lazy <- rvt_data_mapterhorn(pt, res = 2)
  local <- rvt_data_mapterhorn(pt, res = 2, download = TRUE, refresh = TRUE)
  expect_identical(mt_info(lazy)$v, mt_info(local)$v)
  expect_equal(rvt_data_mapterhorn(pt, res = 2, download = TRUE), local)
  whole <- mt_info(rvt_svf(lazy, reach = 20))$v
  pieces <- mt_info(rvt_svf(lazy, reach = 20, tile_size = 137))$v
  expect_identical(whole, pieces)
  # and at native resolution, where the kernel-scale fix alone was not enough
  native <- rvt_data_mapterhorn(pt)
  expect_identical(mt_info(rvt_slope(native))$v,
                   mt_info(rvt_slope(native, tile_size = 137))$v)
})
