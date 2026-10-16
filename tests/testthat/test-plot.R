# rvt_plot()'s stretch. `range` and `pct` are rvt_image()'s spelling, so the
# two functions must resolve the same numbers from the same arguments. Nothing
# here draws: plot_raster() is mocked and the arguments it receives inspected.

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
  svf <- fx_svf()
  a <- capture_plot_args(rvt_plot(svf, pct = c(2, 98)))
  expect_equal(a$minmax_def, rvt_range(svf, pct = c(2, 98)))
})

test_that("a tanh stretch spans the full data and bends the colours, as rvt_image does", {
  m <- fx_msrm()
  full <- rvt_range(m)
  q <- rvtr:::.tanh_range(rvt_range(m, pct = c(2, 98)))   # signed: centred on 0
  a <- capture_plot_args(rvt_plot(m, pct = c(2, 98), stretch = "tanh"))
  # the whole extent is drawn, so nothing is clipped
  expect_equal(a$minmax_def, full)
  # the colour at each data value is the one rvt_image() would give it
  x <- seq(full[1], full[2], length.out = 50)
  got <- grDevices::col2rgb(a$col_map_fn((x - full[1]) / (full[2] - full[1])))[1, ]
  tab <- rvtr:::.col_table("Grays")
  s <- rvtr:::.norm_tanh(x, q[1], q[2], 1.5)
  want <- round(tab[1, pmin(floor(s * 256), 255) + 1] * 255)
  expect_lt(max(abs(got - want)), 3)
  expect_true(is.na(a$col_map_fn(NA_real_)))

  # RGB: one curve across the three bands
  rgb <- fs::file_temp(ext = "tif"); on.exit(unlink(rgb), add = TRUE)
  gdalraster::translate(m, rgb, cl_arg = c("-b", "1", "-b", "1", "-b", "1"), quiet = TRUE)
  a <- capture_plot_args(rvt_plot(rgb, stretch = "tanh"))
  expect_equal(length(a$minmax_def), 6L)
  expect_true(is.function(a$col_map_fn))

  expect_error(rvt_plot(m, stretch = "tanh", minmax_def = c(0, 1)), "not both")
})

test_that("auto picks tanh for signed data and linear for everything else", {
  expect_true(rvtr:::.is_signed(c(-0.21, 0.16)))
  expect_false(rvtr:::.is_signed(c(0.83, 0.99)))           # positive
  expect_false(rvtr:::.is_signed(c(-5, 50)))               # a DEM dipping below sea level
  expect_false(rvtr:::.is_signed(c(0, 1)))                 # touching zero is not crossing it

  tanh_used <- function(p, ...) is.function(capture_plot_args(rvt_plot(p, ...))$col_map_fn)
  expect_true(tanh_used(fx_msrm()))
  expect_false(tanh_used(fx_svf()))
  expect_false(tanh_used(dem_path()))
  expect_false(tanh_used(fx_msrm(), stretch = "linear"))   # naming it overrides
  # signed, but drawn with gdalraster's own stretch: auto steps aside, no error
  a <- capture_plot_args(rvt_plot(fx_msrm(), minmax_pct_cut = c(2, 98)))
  expect_false(is.function(a$col_map_fn))
  expect_equal(a$minmax_pct_cut, c(2, 98))
})

test_that("a signed tanh stretch puts zero on the palette's centre", {
  expect_equal(rvtr:::.tanh_range(c(-0.2, 0.1)), c(-0.2, 0.2))
  expect_equal(rvtr:::.tanh_range(c(0.8, 1)), c(0.8, 1))   # not signed: unchanged
  m <- fx_msrm()
  a <- capture_plot_args(rvt_plot(m, "Blue-Red 3"))
  full <- a$minmax_def
  mid <- a$col_map_fn((0 - full[1]) / (full[2] - full[1]))
  want <- grDevices::colorRamp(grDevices::hcl.colors(256, "Blue-Red 3"))(0.5)
  expect_lt(max(abs(grDevices::col2rgb(mid)[, 1] - want[1, ])), 2)
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
