# Overviews: majority for class codes, cubic for everything else. Read at
# exactly half size, which GDAL serves from the first overview level.

overview_values <- function(p) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  expect_gt(ds$getOverviewCount(1), 0)
  nx <- ds$getRasterXSize() %/% 2; ny <- ds$getRasterYSize() %/% 2
  v <- as.numeric(ds$read(1, 0, 0, nx * 2, ny * 2, nx, ny))
  nd <- ds$getNoDataValue(1)
  if (!is.na(nd)) v[v == nd] <- NA
  v[!is.na(v)]
}

test_that("geomorphons overviews hold only whole class codes", {
  dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
  g <- rvt_geomorphons(dem, reach = 20)
  on.exit(unlink(g), add = TRUE)
  v <- overview_values(g)
  # averaging would invent in-between values that are not classes at all
  expect_true(all(v == round(v)))
  expect_true(all(v %in% rvt_geomorphon_classes))
})

test_that("continuous overviews stay within the data's range beside NoData", {
  dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
  s <- rvt_slope(dem)
  on.exit(unlink(s), add = TRUE)
  v <- overview_values(s)
  # cubic can overshoot, but holes flagged as NoData must not drag
  # neighbouring overview cells towards -9999
  expect_gte(min(v), 0)
  expect_lte(max(v), 90)
})
