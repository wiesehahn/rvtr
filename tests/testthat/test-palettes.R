# Palette names for rvt_plot(), and the palette browser. Everything draws to a
# throwaway PNG so no graphics window opens during the suite.

draw_to_png <- function(expr) {
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  force(expr)
}

test_that("a palette name resolves to that palette's colours", {
  expect_identical(rvtr:::.resolve_col("Grays"), grDevices::hcl.colors(256, "Grays"))
  # case, spaces and hyphens ignored, as hcl.colors() does
  want <- grDevices::hcl.colors(256, "Blue-Red 3")
  expect_identical(rvtr:::.resolve_col("blue-red 3"), want)
  expect_identical(rvtr:::.resolve_col("bluered3"), want)
  # a unique start is enough
  expect_identical(rvtr:::.resolve_col("viri"), grDevices::hcl.colors(256, "Viridis"))
})

test_that("unclear or unknown names fail in plain English", {
  # "red" starts several palettes; hcl.colors() would reject it without saying
  # which, and in the session language
  expect_error(rvtr:::.resolve_col("red"), "matches several palettes")
  expect_error(rvtr:::.resolve_col("notapalette"), "rvt_palettes")
  expect_error(rvtr:::.resolve_col(42), "palette")
})

test_that("a vector of colours is used as given", {
  cols <- c("#000000", "#808080", "#FFFFFF")
  expect_identical(rvtr:::.resolve_col(cols), cols)
  expect_error(rvtr:::.resolve_col(c("black", "notacolour")), "not a colour")
})

test_that("rvt_palettes draws every palette and returns their names", {
  all <- draw_to_png(rvt_palettes())
  expect_identical(all, unlist(lapply(rvtr:::.palette_types, grDevices::hcl.pals),
                               use.names = FALSE))
  div <- draw_to_png(rvt_palettes("diverging"))
  expect_identical(div, grDevices::hcl.pals("diverging"))
  expect_error(rvt_palettes("pastel"), "must be NULL or one of")
})

test_that("rvt_plot takes a palette name as its second argument", {
  dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
  expect_identical(draw_to_png(rvt_plot(dem, "Viridis")), rvtr:::.as_path(dem))
  expect_error(draw_to_png(rvt_plot(dem, "notapalette")), "rvt_palettes")
})
