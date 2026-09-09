# Regression / sanity checks against the bundled rvt-py output, plus the
# structural invariants (mosaic, tiling) that matter regardless of any
# reference implementation.
#
# This package does not aim for bit-for-bit parity with rvt-py - two
# known rvt-py quirks are deliberately not reproduced by default (see
# ?rvt_vat and src/horizon.cpp), and future changes may improve on rvt-py's
# numbers elsewhere too. These tests catch real regressions (wrong sign,
# broken formula, NA leakage) without pinning behaviour to rvt-py's bugs.
# `rvt_compat = TRUE` on rvt_vat() is kept around for anyone who wants to
# verify against rvt-py directly.
#
# The rvt-py comparison outputs (inst/extdata/rvt_*.tif) are gitignored -
# local development aids, not shipped fixtures (see inst/extdata/README.md) -
# so every test that needs one skips itself when the file isn't present.
# dtm1.tif is always committed, so the structural-invariant tests below never
# skip.

ref_path <- function(f) system.file("extdata", f, package = "rvtr")
has_ref <- function(f) nzchar(ref_path(f))

read_band <- function(p) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  m <- matrix(as.numeric(ds$read(1L, 0, 0, nx, ny, nx, ny)),
              nrow = ny, ncol = nx, byrow = TRUE)
  nd <- ds$getNoDataValue(1L)
  if (!is.na(nd)) m[m == nd] <- NA
  m
}

dem <- ref_path("dtm1.tif")

test_that("sky-view factor agrees with rvt-py", {
  skip_if_not(has_ref("rvt_sky_view_factor.tif"), "rvt-py reference not present locally")
  got <- read_band(rvt_svf(dem))
  ref <- read_band(ref_path("rvt_sky_view_factor.tif"))
  expect_equal(dim(got), dim(ref))
  expect_true(all(got >= 0 & got <= 1, na.rm = TRUE))
  expect_gt(cor(as.vector(got), as.vector(ref), use = "complete.obs"), 0.999)
})

test_that("positive openness agrees with rvt-py away from NoData edges", {
  skip_if_not(has_ref("rvt_openness_positive.tif"), "rvt-py reference not present locally")
  got <- read_band(rvt_openness(dem))
  ref <- read_band(ref_path("rvt_openness_positive.tif"))
  expect_equal(dim(got), dim(ref))
  # not bounded to [0, 90]: 90 is flat, above it is convex/elevated terrain
  expect_true(all(got >= 0 & got <= 180, na.rm = TRUE))
  # a handful of pixels beside NoData holes intentionally differ: rvt-py
  # treats "no horizon found" as a slope of -1000 rather than -Inf, off by
  # a fraction of a degree; this package uses the exact value instead.
  expect_gt(cor(as.vector(got), as.vector(ref), use = "complete.obs"), 0.999)
  expect_lt(sum(abs(got - ref) > 0.02, na.rm = TRUE), 5)
})

test_that("rvt_vat(rvt_compat = TRUE) reproduces rvt-py's output", {
  skip_if_not(has_ref("rvt_archaeological_vat_combined.tif"), "rvt-py reference not present locally")
  got <- read_band(rvt_vat(dem, rvt_compat = TRUE))
  ref <- read_band(ref_path("rvt_archaeological_vat_combined.tif"))
  expect_equal(dim(got), dim(ref))
  expect_lt(max(abs(got - ref), na.rm = TRUE), 1e-03)
})

test_that("rvt_vat() default output is sane and strongly correlated with rvt-py", {
  # the default applies the openness layer's declared 50% opacity, which
  # rvt-py's blend_overlay() bug skips - so this is deliberately not a
  # tight match, just a check that nothing is broken.
  skip_if_not(has_ref("rvt_archaeological_vat_combined.tif"), "rvt-py reference not present locally")
  got <- read_band(rvt_vat(dem))
  ref <- read_band(ref_path("rvt_archaeological_vat_combined.tif"))
  expect_equal(dim(got), dim(ref))
  expect_true(all(got >= 0 & got <= 1, na.rm = TRUE))
  expect_gt(cor(as.vector(got), as.vector(ref), use = "complete.obs"), 0.95)
})

test_that("a mosaic gives bit-identical results to the whole raster", {
  # cut the tile into quadrants, then process them as a mosaic
  info <- rvtr:::.dem_info(dem)
  hx <- info$nx %/% 2; hy <- info$ny %/% 2
  quads <- character(0)
  for (yy in list(c(0, hy), c(hy, info$ny - hy))) {
    for (xx in list(c(0, hx), c(hx, info$nx - hx))) {
      q <- tempfile(fileext = ".tif")
      gdalraster::translate(dem, q,
        cl_arg = c("-srcwin", xx[1], yy[1], xx[2], yy[2]), quiet = TRUE)
      quads <- c(quads, q)
    }
  }
  on.exit(unlink(quads), add = TRUE)

  expect_identical(read_band(rvt_svf(quads)), read_band(rvt_svf(dem)))
})

test_that("tiling does not change the result", {
  # a tile size that forces several uneven tiles plus overlap handling
  expect_equal(read_band(rvt_svf(dem, tile_size = 137)),
               read_band(rvt_svf(dem)),
               tolerance = 1e-12)
})
