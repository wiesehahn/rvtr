# rvt_plot()'s stretch. `range` and `pct` are rvt_image()'s spelling, so the
# two functions must resolve the same numbers from the same arguments. Nothing
# here draws: plot_raster() is mocked and the arguments it receives inspected.

dem_path <- function() system.file("extdata", "dtm1.tif", package = "rvtr")

capture_plot_args <- function(expr) {
  args <- NULL
  testthat::with_mocked_bindings(
    force(expr),
    plot_raster = function(...) { args <<- list(...); invisible(NULL) },
    .package = "gdalraster")
  args
}

test_that("range becomes a per-band minmax_def", {
  dem <- dem_path()
  a <- capture_plot_args(rvt_plot(dem, range = c(300, 400)))
  expect_equal(a$minmax_def, c(300, 400))

  # three bands share one range - a stretch per band would shift an RGB
  # image's colour balance
  rgb <- fs::file_temp(ext = "tif")
  on.exit(unlink(rgb), add = TRUE)
  gdalraster::translate(dem, rgb, cl_arg = c("-b", "1", "-b", "1", "-b", "1"),
                        quiet = TRUE)
  a <- capture_plot_args(rvt_plot(rgb, range = c(0, 255)))
  expect_equal(a$nbands, 3L)
  expect_equal(a$minmax_def, c(0, 0, 0, 255, 255, 255))
})

test_that("pct is measured exactly as rvt_image() measures it", {
  svf <- rvt_svf(dem_path())
  on.exit(unlink(svf), add = TRUE)
  a <- capture_plot_args(rvt_plot(svf, pct = c(2, 98)))
  expect_equal(a$minmax_def, rvt_range(svf, pct = c(2, 98)))
})

test_that("conflicting stretch arguments are refused", {
  dem <- dem_path()
  expect_error(rvt_plot(dem, range = c(1, 0)), "lo < hi")
  expect_error(rvt_plot(dem, range = c(0, 1), pct = c(2, 98)), "not both")
  expect_error(rvt_plot(dem, range = c(0, 1), minmax_def = c(0, 1)), "not both")
  # gdalraster's own arguments still work on their own
  a <- capture_plot_args(rvt_plot(dem, minmax_pct_cut = c(2, 98)))
  expect_equal(a$minmax_pct_cut, c(2, 98))
})
