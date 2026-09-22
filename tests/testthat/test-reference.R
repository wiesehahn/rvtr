# Regression / sanity checks against the bundled rvt-py output.
#
# This package does not aim for bit-for-bit parity with rvt-py - several known
# rvt-py quirks are deliberately not reproduced by default (see ?rvt_vat,
# src/horizon.cpp and CLAUDE.md), and future changes may improve on rvt-py
# numbers elsewhere too. These tests catch real regressions (wrong sign,
# broken formula, NA leakage) without pinning behaviour to rvt-py bugs.
# `rvt_compat = TRUE` on rvt_vat() is kept around for anyone who wants to
# verify against rvt-py directly, and is the one thing checked tightly.
#
# The rvt-py comparison outputs (inst/extdata/rvt_*.tif) are gitignored -
# local development aids, not shipped fixtures (see inst/extdata/README.md) -
# so every test that needs one skips itself when the file is not present.
#
# The structural invariants that need no reference implementation moved out
# to test-tiling.R, test-pyramid.R, test-surface.R, test-sun.R and
# test-preprocess.R; shared helpers are in helper-rvtr.R.

test_that("sky-view factor agrees with rvt-py", {
  skip_if_not(has_ref("rvt_sky_view_factor.tif"), "rvt-py reference not present locally")
  got <- read_band(fx_svf())
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

test_that("negative openness agrees with rvt-py away from NoData edges", {
  skip_if_not(has_ref("rvt_openness_negative.tif"), "rvt-py reference not present locally")
  got <- read_band(rvt_openness_negative(dem))
  ref <- read_band(ref_path("rvt_openness_negative.tif"))
  expect_equal(dim(got), dim(ref))
  expect_true(all(got >= 0 & got <= 180, na.rm = TRUE))
  # same -1000-vs--Inf sentinel deviation as positive openness (see above)
  expect_gt(cor(as.vector(got), as.vector(ref), use = "complete.obs"), 0.999)
  expect_lt(sum(abs(got - ref) > 0.02, na.rm = TRUE), 5)
})

test_that("anisotropic sky-view factor agrees with rvt-py", {
  skip_if_not(has_ref("rvt_asvf.tif"), "rvt-py reference not present locally")
  got <- read_band(rvt_asvf(dem))
  # The reference was generated asking rvt-py for asvf_dir = 45, not 315:
  # rvt-py attaches its direction weights to the mirrored azimuth (see the
  # direction test below and CLAUDE.md). Comparing our 315 against its 45
  # therefore checks the weighting maths while leaving the deliberate
  # direction fix out of it - the two agree to float32 rounding.
  ref <- read_band(ref_path("rvt_asvf.tif"))
  expect_equal(dim(got), dim(ref))
  expect_true(all(got >= 0 & got <= 1, na.rm = TRUE))
  expect_lt(max(abs(got - ref), na.rm = TRUE), 1e-05)
})

test_that("anisotropic sky-view factor brightens the direction it is asked to", {
  # No rvt-py reference needed: this checks the property directly, on flat
  # ground with a tall wall filling the eastern third. Cells just west of the
  # wall have their eastern horizon blocked and nothing else, so weighting the
  # eastern sky must cost them the most. North/south would not discriminate -
  # the failure mode this guards against mirrors east and west.
  n <- 120L
  m <- matrix(0, n, n)
  m[, 80:n] <- 50
  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = n,
                           ysize = n, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(0, 1, 0, 0, 0, -1))
  ds$write(1L, 0, 0, n, n, as.vector(t(m)))
  ds$close()
  on.exit(unlink(p), add = TRUE)

  band <- function(x) mean(read_band(x)[40:80, 70:79])
  east  <- band(rvt_asvf(p, main_direction = 90))
  west  <- band(rvt_asvf(p, main_direction = 270))
  north <- band(rvt_asvf(p, main_direction = 0))

  expect_lt(east, north)   # blocked side weighted most -> least sky
  expect_lt(north, west)   # open side weighted most -> most sky
})

test_that("local dominance agrees with rvt-py", {
  skip_if_not(has_ref("rvt_local_dominance.tif"), "rvt-py reference not present locally")
  got <- read_band(rvt_local_dominance(dem))
  ref <- read_band(ref_path("rvt_local_dominance.tif"))
  expect_equal(dim(got), dim(ref))
  # rvt-py accumulates ~264 offset contributions in float32 in place; this
  # package accumulates in double, so a small amount of numerical drift vs.
  # rvt-py is expected (unbiased - verified against a ratio/residual check
  # during development, not just capped here) rather than a sign of a bug.
  expect_gt(cor(as.vector(got), as.vector(ref), use = "complete.obs"), 0.999)
  expect_lt(mean(abs(got - ref), na.rm = TRUE), 0.01)
})

test_that("MSRM agrees with rvt-py", {
  skip_if_not(has_ref("rvt_msrm.tif"), "rvt-py reference not present locally")
  got <- read_band(fx_msrm())
  ref <- read_band(ref_path("rvt_msrm.tif"))
  expect_equal(dim(got), dim(ref))
  expect_gt(cor(as.vector(got), as.vector(ref), use = "complete.obs"), 0.999)
  expect_lt(max(abs(got - ref), na.rm = TRUE), 1e-04)
})

test_that("SLRM agrees with rvt-py", {
  skip_if_not(has_ref("rvt_slrm.tif"), "rvt-py reference not present locally")
  got <- read_band(rvt_slrm(dem))
  ref <- read_band(ref_path("rvt_slrm.tif"))
  expect_equal(dim(got), dim(ref))
  expect_lt(max(abs(got - ref), na.rm = TRUE), 1e-04)
})

test_that("SLRM equals dem minus a plain box mean", {
  # rvt_slrm() takes a shortcut: it reuses the MSRM kernel with radii
  # c(0, radius), relying on a radius-0 box mean being the cell itself. This
  # checks that shortcut against a naive box mean written out longhand, so a
  # change to the integral-image kernel can't quietly redefine SLRM.
  set.seed(1)
  n <- 60L
  m <- matrix(round(runif(n * n, 100, 140), 3), n, n)
  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = n,
                           ysize = n, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(0, 1, 0, 0, 0, -1))
  ds$write(1L, 0, 0, n, n, as.vector(t(m)))
  ds$close()
  on.exit(unlink(p), add = TRUE)

  r <- 5L
  # clamping the window indices replicates the edge cell, which is the same
  # padding the kernel applies
  want <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      ii <- pmin(pmax((i - r):(i + r), 1L), n)
      jj <- pmin(pmax((j - r):(j + r), 1L), n)
      want[i, j] <- m[i, j] - mean(m[ii, jj])
    }
  }

  expect_equal(read_band(rvt_slrm(p, radius = r)), want, tolerance = 1e-5)
})

test_that("slope and hillshade agree with rvt-py", {
  skip_if_not(has_ref("rvt_slope.tif"), "rvt-py reference not present locally")
  got <- read_band(fx_slope())
  ref <- read_band(ref_path("rvt_slope.tif"))
  expect_equal(dim(got), dim(ref))
  expect_true(all(got >= 0 & got <= 90, na.rm = TRUE))
  expect_lt(max(abs(got - ref), na.rm = TRUE), 1e-04)

  got <- read_band(fx_hillshade())
  ref <- read_band(ref_path("rvt_hillshade.tif"))
  expect_true(all(got >= 0 & got <= 1, na.rm = TRUE))
  expect_lt(max(abs(got - ref), na.rm = TRUE), 1e-05)
  # a valid cell next to a NoData hole must still get a value - only the
  # holes themselves are NoData, no eroded ring around them
  expect_equal(sum(is.na(got)), sum(is.na(ref)))
})

test_that("multi-hillshade is one 0-1 band with NoData preserved", {
  p <- rvt_multi_hillshade(dem)
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  n <- ds$getRasterCount()
  ds$close()
  expect_equal(n, 1L)

  m <- read_band(p)
  d <- read_band(dem)
  # GDAL's Byte 1-255 rescaled: strictly above 0, never above 1
  expect_gt(min(m, na.rm = TRUE), 0)
  expect_lte(max(m, na.rm = TRUE), 1)
  # the DEM's holes stay holes, and nothing else is dropped
  expect_identical(is.na(m), is.na(d))

  # lit from all sides, so it must still track a single hillshade closely
  h <- read_band(fx_hillshade())
  ok <- !is.na(m) & !is.na(h)
  expect_gt(stats::cor(m[ok], h[ok]), 0.8)
})

test_that("MSTP agrees with rvt-py", {
  skip_if_not(has_ref("rvt_mstp.tif"), "rvt-py reference not present locally")
  p <- rvt_mstp(dem, local = c(3, 21, 2), meso = c(23, 103, 18),
                broad = c(123, 223, 50))
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  n <- ds$getRasterCount(); ds$close()
  expect_equal(n, 3L)

  for (b in 1:3) {
    got <- read_band(p, b)
    ref <- read_band(ref_path("rvt_mstp.tif"), b)
    expect_true(all(got >= 0 & got <= 1, na.rm = TRUE))
    expect_gt(cor(as.vector(got), as.vector(ref), use = "complete.obs"), 0.999)
  }
  # The finest band divides by a standard deviation that rvt-py computes from
  # float32 squares - catastrophic cancellation at 300 m elevation with a few
  # cm of relief - so a handful of flat-ground pixels differ substantially.
  # Ours is verified against a two-pass mean/sd below; only the bulk is pinned.
  expect_lt(mean(abs(read_band(p, 3L) -
                     read_band(ref_path("rvt_mstp.tif"), 3L)), na.rm = TRUE), 0.01)
})

test_that("MSTP's deviation survives large elevation offsets", {
  # The summed-area variance is the numerically delicate part: at 3000 m base
  # elevation with centimetre relief, E[z^2] and E[z]^2 agree to many digits
  # before differing. Checked against a two-pass mean/sd, which is stable.
  n <- 40L
  set.seed(7)
  relief <- matrix(round(runif(n * n, -0.05, 0.05), 4), n, n)
  pad <- 6L
  for (base in c(0, 3000)) {
    padded <- rvtr:::.pad_symmetric(base + relief, pad)
    got <- rvtr:::max_deviation_kernel(padded, pad, n, n, 3L, 3L, 1L, 1L)
    want <- matrix(NA_real_, n, n)
    for (i in seq_len(n)) for (j in seq_len(n)) {
      w <- as.vector(padded[(i - 3):(i + 3) + pad, (j - 3):(j + 3) + pad])
      mu <- mean(w)
      want[i, j] <- (padded[i + pad, j + pad] - mu) /
        (sqrt(mean((w - mu)^2)) + 1e-6)
    }
    expect_equal(got, want, tolerance = 1e-8)
  }
})

test_that("sky illumination and shadow read 1 on flat open ground", {
  # No reference file needed. Flat ground has no horizon anywhere, so every
  # direction contributes its full share and the flat-ground normalisation
  # must land exactly on 1; and nothing can shadow it.
  n <- 80L
  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = n,
                           ysize = n, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(0, 1, 0, 0, 0, -1))
  ds$write(1L, 0, 0, n, n, rep(100, n * n))
  ds$close()
  on.exit(unlink(p), add = TRUE)

  for (model in c("uniform", "overcast")) {
    v <- read_band(rvt_sky_illumination(p, sky_model = model, reach = 8,
                                        num_directions = 16, overwrite = TRUE))
    expect_equal(range(v), c(1, 1), tolerance = 1e-9,
                 info = paste("sky model", model))
  }
  expect_true(all(read_band(rvt_shadow(p, sun_elevation = 15, reach = 8)) == 1))
})

test_that("cast shadow agrees with rvt-py", {
  skip_if_not(has_ref("rvt_shadow.tif"), "rvt-py reference not present locally")
  skip_if_not(has_ref("rvt_simfilled_dem.tif"), "rvt-py reference not present locally")
  # rvt-py's sky_illumination() cannot run on a DEM containing NoData at all
  # (it pads with np.min, which is NaN), so the reference uses a hole-filled
  # copy; compare on the same input.
  filled <- ref_path("rvt_simfilled_dem.tif")
  got <- read_band(rvt_shadow(filled, sun_azimuth = 315, sun_elevation = 15,
                              reach = 30))
  ref <- read_band(ref_path("rvt_shadow.tif"))
  expect_true(all(got %in% c(0, 1)))
  # rvt-py snaps the sun to its nearest search direction and approximates the
  # horizon through a pyramid; we trace the exact azimuth at full resolution
  expect_gt(mean(got == ref), 0.99)
})

