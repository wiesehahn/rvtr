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

## A cropped copy, for the tests that compare the multi-resolution search
## against an exact full-resolution one. That baseline - pyramid_px equal to
## the whole reach - is what costs, not the pyramid, so the saving comes from
## the raster rather than from weakening the approximation: pyramid_px has to
## stay at its default or the error being measured is a different one.
## The interior margin used with these has to *exceed* the reach, not equal
## it, or the comparison picks up the edge extrapolation the two methods
## deliberately differ on: at 600 px with a 200 px margin the sky-view error
## reads 1.2e-4 against 3.4e-5 well clear of the edge.
crop_raster <- function(src, size = 800L) {
  p <- tempfile(fileext = ".tif")
  gdalraster::translate(src, p, quiet = TRUE,
                         cl_arg = c("-srcwin", "0", "0", as.character(size),
                                     as.character(size)))
  p
}

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
  got <- read_band(rvt_msrm(dem))
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
  got <- read_band(rvt_slope(dem))
  ref <- read_band(ref_path("rvt_slope.tif"))
  expect_equal(dim(got), dim(ref))
  expect_true(all(got >= 0 & got <= 90, na.rm = TRUE))
  expect_lt(max(abs(got - ref), na.rm = TRUE), 1e-04)

  got <- read_band(rvt_hillshade(dem))
  ref <- read_band(ref_path("rvt_hillshade.tif"))
  expect_true(all(got >= 0 & got <= 1, na.rm = TRUE))
  expect_lt(max(abs(got - ref), na.rm = TRUE), 1e-05)
  # a valid cell next to a NoData hole must still get a value - only the
  # holes themselves are NoData, no eroded ring around them
  expect_equal(sum(is.na(got)), sum(is.na(ref)))
})

read_band_n <- function(p, b) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  m <- matrix(as.numeric(ds$read(b, 0, 0, nx, ny, nx, ny)),
              nrow = ny, ncol = nx, byrow = TRUE)
  nd <- ds$getNoDataValue(b)
  if (!is.na(nd)) m[m == nd] <- NA
  m
}

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
  h <- read_band(rvt_hillshade(dem))
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
    got <- read_band_n(p, b)
    ref <- read_band_n(ref_path("rvt_mstp.tif"), b)
    expect_true(all(got >= 0 & got <= 1, na.rm = TRUE))
    expect_gt(cor(as.vector(got), as.vector(ref), use = "complete.obs"), 0.999)
  }
  # The finest band divides by a standard deviation that rvt-py computes from
  # float32 squares - catastrophic cancellation at 300 m elevation with a few
  # cm of relief - so a handful of flat-ground pixels differ substantially.
  # Ours is verified against a two-pass mean/sd below; only the bulk is pinned.
  expect_lt(mean(abs(read_band_n(p, 3L) -
                     read_band_n(ref_path("rvt_mstp.tif"), 3L)), na.rm = TRUE), 0.01)
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

make_raster <- function(m) {
  n_r <- nrow(m); n_c <- ncol(m)
  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = n_c,
                           ysize = n_r, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(0, 1, 0, 0, 0, -1))
  ds$write(1L, 0, 0, n_c, n_r, as.vector(t(m)))
  ds$close()
  p
}

test_that("aspect points down the slope", {
  n <- 120L
  # plane falling to the east -> aspect 90; falling to the south -> 180
  east  <- make_raster(matrix(100 - 0.1 * col(matrix(0, n, n)), n, n))
  south <- make_raster(matrix(100 - 0.1 * row(matrix(0, n, n)), n, n))
  on.exit(unlink(c(east, south)), add = TRUE)
  mid <- function(x) median(read_band(x)[20:100, 20:100])
  expect_equal(mid(rvt_aspect(east)), 90, tolerance = 1e-6)
  expect_equal(mid(rvt_aspect(south)), 180, tolerance = 1e-6)
})

test_that("rvt_tpi is rvt_slrm under its other name", {
  expect_identical(read_band(rvt_tpi(dem, radius = 7)),
                   read_band(rvt_slrm(dem, radius = 7)))
})

test_that("DEV is SLRM standardised by local roughness", {
  # flat ground has no deviation and no roughness; the 1e-6 guard keeps it
  # finite rather than 0/0
  flat <- make_raster(matrix(100, 60L, 60L))
  on.exit(unlink(flat), add = TRUE)
  expect_true(all(abs(read_band(rvt_dev(flat, radius = 5))) < 1e-6))

  r <- 10L
  d <- read_band(rvt_dev(dem, radius = r))
  expect_true(all(is.finite(d[!is.na(d)])))

  # Away from the border DEV must agree in sign with SLRM exactly - it is the
  # same numerator over a positive divisor. Within `radius` of the edge the
  # two legitimately differ: DEV mirrors the terrain outwards (so the
  # standard deviation it divides by doesn't collapse) while SLRM replicates
  # the edge, giving different window means.
  s <- read_band(rvt_slrm(dem, radius = r))
  n <- nrow(d)
  edge <- pmin(pmin(row(d), n - row(d) + 1L), pmin(col(d), n - col(d) + 1L))
  ok <- !is.na(d) & !is.na(s) & abs(s) > 1e-6 & edge > r
  expect_identical(sign(d[ok]), sign(s[ok]))
})

test_that("geomorphons classify the obvious shapes", {
  n <- 101L; c0 <- 51L
  rad <- sqrt((row(matrix(0, n, n)) - c0)^2 + (col(matrix(0, n, n)) - c0)^2)
  cone <- make_raster(matrix(100 + pmax(0, 20 - 0.5 * rad), n, n))
  pit  <- make_raster(matrix(100 - pmax(0, 20 - 0.5 * rad), n, n))
  flat <- make_raster(matrix(100, n, n))
  on.exit(unlink(c(cone, pit, flat)), add = TRUE)

  cls <- rvt_geomorphon_classes
  expect_equal(read_band(rvt_geomorphons(cone, reach = 15))[c0, c0],
               as.numeric(cls[["peak"]]))
  expect_equal(read_band(rvt_geomorphons(pit, reach = 15))[c0, c0],
               as.numeric(cls[["pit"]]))
  # perfectly level ground: every direction is within the flatness threshold
  f <- read_band(rvt_geomorphons(flat, reach = 15))
  expect_true(all(f[30:70, 30:70] == as.numeric(cls[["flat"]])))

  # real terrain: only valid class codes, and NoData preserved
  g <- read_band(rvt_geomorphons(dem, reach = 10))
  expect_true(all(g[!is.na(g)] %in% as.numeric(cls)))
  expect_equal(sum(is.na(g)), sum(is.na(read_band(dem))))
})

test_that("curvature matches analytic surfaces", {
  n <- 81L; mid <- 30:52
  x <- col(matrix(0, n, n)) - 41
  y <- row(matrix(0, n, n)) - 41
  cur <- function(p, ty) read_band(rvt_curvature(p, type = ty))[mid, mid]

  # a plane has no curvature of any kind
  pl <- make_raster(matrix(100 + 0.2 * x + 0.1 * y, n, n))
  for (ty in c("profile", "plan", "mean", "total", "gaussian"))
    expect_lt(max(abs(cur(pl, ty))), 1e-4, label = ty)

  # paraboloid z = a(x^2 + y^2): zxx = zyy = 2a, so mean = 4a, gaussian = 4a^2.
  # These pin the scale - the Zevenbergen & Thorne polynomial coefficients are
  # half the partial derivatives, which is easy to conflate.
  # A *tilted* paraboloid is the load-bearing case: with a non-zero gradient
  # the proper curvatures diverge from the bare second-derivative
  # combinations, so this catches a missing (1+slope^2) denominator. An
  # untilted one cannot - the two agree there, which is how that bug survived
  # its first test.
  a <- 0.01; tilt <- 0.3
  pb <- make_raster(matrix(100 + a * (x^2 + y^2) + tilt * x, n, n))
  fxx <- 2 * a; fyy <- 2 * a; q <- 1 + tilt^2
  H <- -((1 + tilt^2) * fyy + fxx) / (2 * q^1.5)   # convex-positive
  K <- (fxx * fyy) / q^2
  d <- sqrt(max(0, H^2 - K))
  inner <- 38:44
  cur2 <- function(p, ty) read_band(rvt_curvature(p, type = ty))[inner, inner]
  expect_equal(median(cur2(pb, "mean")), H, tolerance = 1e-2)
  expect_equal(median(cur2(pb, "gaussian")), K, tolerance = 1e-2)
  expect_equal(median(cur2(pb, "minimal")), H - d, tolerance = 1e-2)
  expect_equal(median(cur2(pb, "maximal")), H + d, tolerance = 1e-2)
  # total follows Evans' un-normalised definition, matching other software
  expect_equal(median(cur2(pb, "total")), fxx^2 + fyy^2, tolerance = 1e-2)

  # saddle z = c*x*y: mean 0, gaussian -c^2. Only this exercises the mixed
  # partial, so it is what catches a sign error there.
  cc <- 0.01
  sd_ <- make_raster(matrix(100 + cc * x * y, n, n))
  expect_lt(abs(median(cur(sd_, "mean"))), 1e-6)
  expect_equal(median(cur(sd_, "gaussian")), -cc^2, tolerance = 1e-2)
  # a saddle has principal curvatures of opposite sign
  expect_lt(median(cur(sd_, "minimal")), 0)
  expect_gt(median(cur(sd_, "maximal")), 0)

  # sign convention: convex positive, concave negative, on a diagonal slope
  # too (where the cross term contributes)
  ridge  <- make_raster(matrix(100 - 0.02 * x^2, n, n))
  trough <- make_raster(matrix(100 + 0.02 * ((x + y) / sqrt(2))^2 + 0.3 * x, n, n))
  expect_gt(median(cur(ridge, "profile")), 0)
  expect_lt(median(cur(trough, "profile")), 0)

  unlink(c(pl, pb, sd_, ridge, trough))
})

test_that("Laplacian-of-Gaussian is flat on a plane and signed across an edge", {
  n <- 81L
  x <- col(matrix(0, n, n)) - 41
  pl <- make_raster(matrix(100 + 0.2 * x, n, n))
  expect_lt(max(abs(read_band(rvt_log(pl, sigma = 2))[30:52, 30:52])), 1e-6)

  step <- matrix(100, n, n); step[, 42:n] <- 101
  sp <- make_raster(step)
  lg <- read_band(rvt_log(sp, sigma = 2))[41, ]
  # a step produces the classic antisymmetric pair: positive on the low side,
  # negative on the high side, summing to ~0
  expect_gt(max(lg[35:41]), 0)
  expect_lt(min(lg[42:48]), 0)
  expect_equal(sum(lg[35:48]), 0, tolerance = 1e-6)

  unlink(c(pl, sp))
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

test_that("tiling does not change the MSRM result", {
  # separate from the horizon-kernel check above: MSRM uses a different
  # kernel (integral image, not horizon search), worth its own invariant
  expect_equal(read_band(rvt_msrm(dem, tile_size = 137)),
               read_band(rvt_msrm(dem)),
               tolerance = 1e-12)
})

test_that("rvt_resample only coarsens, and coarsens correctly", {
  p <- rvt_resample(dem, 4)
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  nx <- ds$getRasterXSize(); r <- ds$res()
  ds$close()
  expect_equal(nx, 250L)
  expect_equal(unname(r[1]), 4)

  # detail cannot be invented
  expect_error(rvt_resample(dem, 0.5), "finer than the source")
  # asking for the source resolution is allowed - a no-op copy, same grid
  same <- rvt_resample(dem, 1)
  ds <- methods::new(gdalraster::GDALRaster, same, read_only = TRUE)
  expect_equal(ds$getRasterXSize(), 1000L)
  ds$close()

  # the elevations must stay elevations - no NoData bled into the mean
  src <- read_band(dem)
  got <- read_band(p)
  expect_gte(min(got, na.rm = TRUE), min(src, na.rm = TRUE))
  expect_lte(max(got, na.rm = TRUE), max(src, na.rm = TRUE))
})

test_that("a metric on a resampled DEM is still tile-invariant", {
  # the reason rvt_resample() writes a real file instead of a lazy VRT:
  # GDAL's on-the-fly resampling can be window-dependent, which would make
  # results depend on tile_size. Both an integer and a non-integer factor.
  for (res in c(4, 2.5)) {
    p <- rvt_resample(dem, res)
    expect_identical(read_band(rvt_svf(p, tile_size = 37)),
                     read_band(rvt_svf(p)))
  }
})

test_that("categorical resampling keeps whole class codes", {
  # pins `-ovr NONE`: without it GDAL serves the source COG's AVERAGE-built
  # overviews and silently ignores method = "mode", returning class codes
  # like 4.31 that correspond to no geomorphon at all.
  g <- rvt_geomorphons(dem)
  classes <- sort(unique(as.vector(read_band(g))))
  classes <- classes[!is.na(classes)]

  m <- read_band(rvt_resample(g, 4, method = "mode"))
  got <- sort(unique(as.vector(m)))
  got <- got[!is.na(got)]
  expect_true(all(got == round(got)))
  expect_true(all(got %in% classes))
})

test_that("solar positions match known geometry", {
  # noon altitude at the solstice is 90 - latitude + declination
  s <- rvtr:::.solar_positions(as.Date("2025-06-21"), 1, 51.44, 0)
  expect_equal(max(s$altitude), 90 - 51.44 + 23.44, tolerance = 0.01)
  # the equator gets ~12 h whatever the season
  for (d in c("2025-03-20", "2025-06-21", "2025-12-21")) {
    s <- rvtr:::.solar_positions(as.Date(d), 1, 0, 0)
    expect_equal(sum(s$altitude > 0) / 60, 12, tolerance = 0.02)
  }
  # at an equinox the sun rises due east
  s <- rvtr:::.solar_positions(as.Date("2025-03-20"), 1, 51.44, 0)
  up <- which(s$altitude > 0)
  expect_equal(s$azimuth[up[1]], 90, tolerance = 0.5)
})

test_that("time in daylight is exact on open flat ground", {
  flat <- make_raster(matrix(100, 60, 60))
  day <- as.Date("2025-06-21")
  lat <- 51.44; lon <- 9.83
  expected <- sum(rvtr:::.solar_positions(day, 10, lat, lon)$altitude > 0) * 10 / 60

  # tolerance is float32 rounding: results are stored as Float32
  h <- read_band(rvt_daylight(flat, dates = day, lat = lat, lon = lon,
                               reach = 20))
  expect_equal(max(h), expected, tolerance = 1e-6)
  expect_equal(min(h), expected, tolerance = 1e-6)

  f <- read_band(rvt_daylight(flat, dates = day, lat = lat, lon = lon,
                               reach = 20, units = "fraction"))
  expect_equal(max(f), 1, tolerance = 1e-6)
})

test_that("a shaft gets no sun while the plateau around it gets all of it", {
  m <- matrix(500, 60, 60)
  m[28:32, 28:32] <- 0
  p <- make_raster(m)
  day <- as.Date("2025-06-21")
  h <- read_band(rvt_daylight(p, dates = day, lat = 51.44, lon = 9.83,
                               reach = 20))
  expect_equal(h[30, 30], 0)
  expect_gt(h[5, 5], 16)
})

test_that("slope aspect and season interact the way the sun does", {
  # A 30 degree slope at 51 N. In midsummer the sun rises in the north-east
  # and is 62 degrees high at noon, so the north-facing slope is lit all day
  # while the south-facing one is self-shadowed at dawn and dusk. Come the
  # winter solstice the sun peaks at 15 degrees and the ordering inverts
  # completely. Counting hours is not the same as counting energy.
  n <- 80
  fall <- tan(30 * pi / 180)
  south <- outer(seq(n, 1) * fall, rep(1, n))   # drops southward
  north <- outer(seq(1, n) * fall, rep(1, n))   # drops northward
  hours <- function(m, d) {
    read_band(rvt_daylight(make_raster(m), dates = as.Date(d), lat = 51.44,
                            lon = 9.83, reach = 20))[40, 40]
  }
  expect_gt(hours(north, "2025-06-21"), hours(south, "2025-06-21"))
  expect_equal(hours(north, "2025-12-21"), 0)
  expect_gt(hours(south, "2025-12-21"), 7)
})

test_that("tiling does not change time in daylight", {
  args <- list(dem, dates = as.Date("2025-06-21"), reach = 20, time_step = 30)
  expect_equal(read_band(do.call(rvt_daylight, c(args, list(tile_size = 137)))),
               read_band(do.call(rvt_daylight, args)),
               tolerance = 1e-12)
})

test_that("the pyramid plan covers the requested reach", {
  p <- .pyramid_plan(xres = 1, reach = 400, px = 100, factor = 4)
  expect_equal(length(p), 2L)
  expect_equal(p[[1]]$fac, 1L)          # level 0 is always the DEM itself
  expect_equal(p[[1]]$outer, 100)
  expect_equal(p[[2]]$fac, 4L)
  expect_equal(p[[2]]$outer, 400)       # reaches exactly as far as asked
  expect_equal(p[[2]]$inner, 100)       # and starts where level 0 stopped

  # levels stay contiguous and the last one always lands on `reach`
  for (cfg in list(c(0.25, 1600, 100, 4), c(1, 250, 50, 2), c(2, 5000, 100, 8))) {
    p <- .pyramid_plan(cfg[1], cfg[2], cfg[3], cfg[4])
    expect_equal(p[[length(p)]]$outer, cfg[2])
    for (i in seq_along(p)[-1]) expect_equal(p[[i]]$inner, p[[i - 1]]$outer)
  }
  expect_error(.pyramid_plan(1, -5, 100, 4), "positive")
  expect_error(.pyramid_plan(1, 400, 100, 1), "at least 2")
})

test_that("a multi-resolution search matches a full-resolution one", {
  # Compared away from the raster edge: within `reach` of the boundary both
  # methods are extrapolating terrain that isn't there, and they do it
  # differently (level 0 reflects, coarse levels replicate).
  d <- crop_raster(dem)
  on.exit(unlink(d), add = TRUE)
  interior <- function(m, k) m[(k + 1):(nrow(m) - k), (k + 1):(ncol(m) - k)]
  err <- function(a, b, k = 300) {
    a <- interior(read_band(a), k); b <- interior(read_band(b), k)
    ok <- !is.na(a) & !is.na(b)
    mean(abs(a[ok] - b[ok]))
  }
  # reach 200 with the default pyramid_px gives two levels; pyramid_px = 200
  # covers the same reach in one, which is the exact search to compare against
  expect_lt(err(rvt_svf(d, reach = 200, pyramid_px = 200),
                rvt_svf(d, reach = 200)), 1e-4)
  expect_lt(err(rvt_openness(d, reach = 200, pyramid_px = 200),
                rvt_openness(d, reach = 200)), 0.02)
  # Negative openness inverts every level, which only works if the coarse
  # levels take the mirror statistic (min for max - see .mirror_method()).
  # Its tolerance is deliberately looser than positive openness's 0.02: on
  # this terrain the error is 0.020 against 0.004, and it does not shrink with
  # a wider margin, so it is not an edge effect. Inverting turns narrow
  # incised gullies into narrow ridges, and narrow features are exactly what
  # coarsening represents least well - the same asymmetry the distant-wall
  # measurements show.
  expect_lt(err(rvt_openness_negative(d, reach = 200, pyramid_px = 200),
                rvt_openness_negative(d, reach = 200)), 0.03)
})

test_that("a coarser pyramid is worse but still close", {
  # Reference and candidates must share a reach - the comparison is about
  # pyramid_px alone. 200 m keeps a far field worth measuring while costing a
  # quarter of what 400 m did.
  d <- crop_raster(dem)
  on.exit(unlink(d), add = TRUE)
  interior <- function(m, k = 300) m[(k + 1):(nrow(m) - k), (k + 1):(ncol(m) - k)]
  ref <- interior(read_band(rvt_openness(d, reach = 200, pyramid_px = 200)))
  e <- function(px) {
    got <- interior(read_band(rvt_openness(d, reach = 200, pyramid_px = px)))
    ok <- !is.na(got) & !is.na(ref)
    mean(abs(got[ok] - ref[ok]))
  }
  # pyramid_px is the accuracy knob: more full-resolution reach, less error
  expect_lt(e(100), e(50))
  expect_lt(e(50), e(25))
})

test_that("tiling and mosaicking do not change a multi-resolution result", {
  one <- read_band(rvt_svf(dem, reach = 100, pyramid_px = 25))
  expect_identical(read_band(rvt_svf(dem, reach = 100, pyramid_px = 25,
                                      tile_size = 137)), one)

  quads <- vapply(1:4, function(i) {
    p <- tempfile(fileext = ".tif")
    gdalraster::translate(dem, p, quiet = TRUE, cl_arg = c("-srcwin",
      as.character(c(0, 500)[1 + (i - 1) %% 2]),
      as.character(c(0, 500)[1 + (i - 1) %/% 2]), "500", "500"))
    p
  }, "")
  on.exit(unlink(quads), add = TRUE)
  expect_identical(read_band(rvt_svf(quads, reach = 100,
                                      pyramid_px = 25)), one)
})

test_that("distances snap to whole cells by a stated rule", {
  # halfway cases must not inherit round()'s round-half-to-even, which would
  # send 2.5 down and 3.5 up
  expect_equal(.cells(5, 2, "d", "outer"), 3L)     # 2.5 -> 3, not 2
  expect_equal(.cells(7, 2, "d", "outer"), 4L)     # 3.5 -> 4
  expect_equal(.cells(4, 2, "d", "outer"), 2L)
  expect_equal(.cells(10, 1, "d", "outer"), 10L)   # exact, and silent
  expect_silent(.cells(10, 1, "d", "outer"))

  # below one cell: an outer radius cannot be honoured at all
  expect_error(.cells(1, 2, "radius", "outer"), "smaller than one cell")
  # but an inner radius or skip may legitimately be zero
  expect_equal(.cells(1, 2, "skip", "inner"), 1L)
  expect_equal(.cells(0.2, 2, "skip", "inner"), 0L)
  # and a step finer than the grid just means every cell
  expect_equal(.cells(0.2, 2, "step", "step"), 1L)

  # say so when the grid could not honour the request closely
  expect_message(.cells(5, 2, "radius", "outer"), "snapped to 3 cells")
})

test_that("metric parameters are distances, not pixel counts", {
  # a 2 m raster: the same request must mean the same ground distance
  n <- 120L
  m <- matrix(rnorm(n * n, 100, 2), n, n)
  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = n,
                           ysize = n, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(0, 2, 0, 0, 0, -2))        # 2 m cells
  ds$write(1L, 0, 0, n, n, as.vector(t(m)))
  ds$close()
  on.exit(unlink(p), add = TRUE)

  # asking for less than one cell is refused rather than silently rounded
  expect_error(rvt_svf(p, reach = 1), "smaller than one cell")
  expect_error(rvt_slrm(p, radius = 1), "smaller than one cell")
  # and a request the grid cannot hit exactly says so
  expect_message(rvt_svf(p, reach = 5), "snapped")

  # 20 m of smoothing on a 2 m raster must equal 10 cells, which is what the
  # same 20 m means on the 1 m sample tile
  expect_equal(read_band(rvt_slrm(p, radius = 20)),
               read_band(rvt_slrm(p, radius = 20.4)), tolerance = 1e-12)
})

test_that("the day range selects the days it says it does", {
  yr <- as.integer(format(Sys.Date(), "%Y"))
  jan1 <- as.Date(paste0(yr, "-01-01"))
  # a stand-in for the date building inside rvt_daylight()
  span <- function(...) {
    f <- formals(rvt_daylight)
    a <- list(...)
    n_year <- as.integer(as.Date(paste0(yr + 1L, "-01-01")) - jan1)
    s <- if (is.null(a$start_day)) f$start_day else a$start_day
    e <- if (is.null(a$end_day)) f$end_day else a$end_day
    st <- if (is.null(a$day_step)) f$day_step else a$day_step
    doy <- if (e >= s) seq(s, e) else c(seq(s, n_year), seq_len(e))
    jan1 + (doy[seq(1, length(doy), by = st)] - 1L)
  }
  expect_equal(length(span()), 73L)                  # whole year, every 5th day
  expect_equal(min(span()), jan1)
  expect_equal(length(span(start_day = 1, end_day = 10, day_step = 1)), 10L)
  # a range ending before it starts wraps through New Year
  w <- span(start_day = 360, end_day = 5, day_step = 1)
  expect_equal(length(w), 11L)
  expect_true(any(format(w, "%m") == "12") && any(format(w, "%m") == "01"))

  expect_error(rvt_daylight(dem, start_day = 0), "at least 1")
  expect_error(rvt_daylight(dem, day_step = 2.5), "whole number")
})

test_that("a narrower day range gives a different answer than the full year", {
  # cheap settings: a handful of days, coarse time step, short reach
  june <- read_band(rvt_daylight(dem, start_day = 152, end_day = 182,
                                  day_step = 15, time_step = 60, reach = 10))
  dec <- read_band(rvt_daylight(dem, start_day = 335, end_day = 365,
                                 day_step = 15, time_step = 60, reach = 10))
  expect_gt(mean(june, na.rm = TRUE), mean(dec, na.rm = TRUE))
  # midsummer must beat midwinter everywhere at this latitude
  expect_true(all(june >= dec, na.rm = TRUE))
})

test_that("daylight, shadow and sky illumination also search multi-resolution", {
  # Each of these has its own kernel, all three now driven by the same
  # RvtLevels accumulation (src/levels.h). Compared against an exact
  # full-resolution search of the same reach - which is what a pyramid_px
  # large enough to cover the whole reach in one level gives - away from the
  # raster edge, where the two extrapolate differently.
  d <- crop_raster(dem)
  on.exit(unlink(d), add = TRUE)
  interior <- function(m, k = 300) m[(k + 1):(nrow(m) - k), (k + 1):(ncol(m) - k)]
  err <- function(a, b) {
    a <- interior(read_band(a)); b <- interior(read_band(b))
    ok <- !is.na(a) & !is.na(b)
    mean(abs(a[ok] - b[ok]))
  }

  expect_lt(err(rvt_shadow(d, reach = 200, pyramid_px = 200),
                rvt_shadow(d, reach = 200)), 0.01)

  expect_lt(err(rvt_sky_illumination(d, num_directions = 8, reach = 200,
                                      pyramid_px = 200),
                rvt_sky_illumination(d, num_directions = 8, reach = 200)), 1e-3)

  expect_lt(err(rvt_daylight(d, dates = as.Date("2025-06-21"), time_step = 60,
                              num_directions = 8, reach = 200, pyramid_px = 200),
                rvt_daylight(d, dates = as.Date("2025-06-21"), time_step = 60,
                              num_directions = 8, reach = 200)), 0.05)
})

test_that("tiling does not change a multi-resolution daylight or shadow", {
  one <- read_band(rvt_shadow(dem, reach = 100, pyramid_px = 25))
  expect_identical(read_band(rvt_shadow(dem, reach = 100, pyramid_px = 25,
                                         tile_size = 137)), one)

  args <- list(dem, dates = as.Date("2025-06-21"), time_step = 60,
               num_directions = 8, reach = 100, pyramid_px = 25)
  expect_identical(read_band(do.call(rvt_daylight, c(args, list(tile_size = 137)))),
                   read_band(do.call(rvt_daylight, args)))

  si <- list(dem, num_directions = 8, reach = 100, pyramid_px = 25)
  expect_identical(read_band(do.call(rvt_sky_illumination, c(si, list(tile_size = 137)))),
                   read_band(do.call(rvt_sky_illumination, si)))
})

# Four adjacent 1 km tiles around the Partnachklamm - 494 m of relief, against
# dtm1.tif's 71 m, and a genuine four-file mosaic rather than one raster cut
# into quadrants. Gitignored like the rvt-py references; see
# inst/extdata/README.md for the one-command download.
alp_tiles <- c("dgm1_659_5259.tif", "dgm1_660_5259.tif",
               "dgm1_659_5258.tif", "dgm1_660_5258.tif")
has_alps <- all(vapply(alp_tiles, has_ref, logical(1)))
alps <- function() vapply(alp_tiles, ref_path, "")

test_that("a mosaic of four separate files matches the merged raster", {
  skip_if_not(has_alps, "alpine tiles not present locally")
  merged <- tempfile(fileext = ".tif")
  gdalraster::translate(rvt_mosaic(alps()), merged, quiet = TRUE)
  on.exit(unlink(merged), add = TRUE)

  # the ordinary single-resolution search
  expect_identical(read_band(rvt_svf(alps())), read_band(rvt_svf(merged)))

  # and the multi-resolution one, whose coarse levels are built from the VRT
  # rather than from a real file - the case quadrant-splitting cannot reach
  expect_identical(read_band(rvt_svf(alps(), reach = 100, pyramid_px = 25)),
                   read_band(rvt_svf(merged, reach = 100, pyramid_px = 25)))
})

test_that("the multi-resolution search holds up in mountains", {
  skip_if_not(has_alps, "alpine tiles not present locally")
  one <- crop_raster(ref_path("dgm1_659_5258.tif"))
  on.exit(unlink(one), add = TRUE)
  interior <- function(m, k = 300) m[(k + 1):(nrow(m) - k), (k + 1):(ncol(m) - k)]

  near <- interior(read_band(rvt_openness(one, reach = 25)))
  far <- interior(read_band(rvt_openness(one, reach = 200, pyramid_px = 200)))
  pyr <- interior(read_band(rvt_openness(one, reach = 200)))

  # the far field has to matter here, or the comparison below proves nothing
  expect_gt(mean(abs(near - far)), 1)

  # and the pyramid must recover nearly all of it
  expect_lt(mean(abs(pyr - far)), 0.1 * mean(abs(near - far)))
})

test_that("daylight in a gorge is far shorter than on the ridge above it", {
  skip_if_not(has_alps, "alpine tiles not present locally")
  one <- ref_path("dgm1_659_5258.tif")
  h <- read_band(rvt_daylight(one, dates = as.Date("2025-12-21"), time_step = 30,
                               num_directions = 16, reach = 500))
  z <- read_band(one)
  # midwinter: the lowest ground in a 494 m gorge must lose hours the tops keep
  low <- h[z < stats::quantile(z, 0.05)]
  high <- h[z > stats::quantile(z, 0.95)]
  expect_lt(mean(low), mean(high))
  expect_lt(min(h), 1)                 # somewhere gets almost no sun at all
  expect_true(all(h >= 0 & h <= 9))    # day length at 47.5 N is about 8.6 h
})

test_that("total hours is the sum over the period, mean hours the daily average", {
  flat <- make_raster(matrix(100, 60, 60))
  lat <- 51.44; lon <- 9.83
  # a short, fully sampled span so the total is unambiguous
  args <- list(flat, start_day = 100, end_day = 109, day_step = 1,
               time_step = 30, lat = lat, lon = lon, reach = 20)

  mean_h <- read_band(do.call(rvt_daylight, c(args, list(units = "mean_hours"))))
  total_h <- read_band(do.call(rvt_daylight, c(args, list(units = "total_hours"))))
  expect_equal(max(total_h), max(mean_h) * 10, tolerance = 1e-5)

  # on open flat ground the total must equal the astronomical day lengths added up
  dates <- as.Date(paste0(format(Sys.Date(), "%Y"), "-01-01")) + (100:109 - 1)
  expected <- sum(rvtr:::.solar_positions(dates, 30, lat, lon)$altitude > 0) * 30 / 60
  expect_equal(max(total_h), expected, tolerance = 1e-5)

  # sampling fewer days must not change what "total" means, only its precision
  sparse <- read_band(do.call(rvt_daylight,
    c(args[names(args) != "day_step"], list(day_step = 5, units = "total_hours"))))
  expect_equal(max(sparse), max(total_h), tolerance = 0.02)
})

test_that("max resampling really takes the block maximum", {
  # gdal_translate cannot do -r max and silently substitutes nearest, which
  # would make the pyramid's coarse levels miss narrow distant obstructions
  # entirely; rvt_resample routes these methods through gdalwarp instead.
  z <- matrix(0, 64, 64)
  z[30:31, ] <- 50                      # a 2-cell ridge, thinner than one block
  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = 64,
                           ysize = 64, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(557000, 1, 0, 5700000, 0, -1))
  ds$setProjection(gdalraster::srs_to_wkt("EPSG:25832"))
  ds$write(1L, 0, 0, 64, 64, as.vector(t(z)))
  ds$close()
  on.exit(unlink(p), add = TRUE)

  got <- read_band(rvt_resample(p, 8, method = "max"))
  expect_equal(max(got), 50)                  # the ridge survives
  expect_equal(sum(got == 50), 8)             # in exactly one row of blocks
  # nearest, by contrast, steps straight over it
  expect_equal(max(read_band(rvt_resample(p, 8, method = "near"))), 0)

  # and min is the mirror image, which is what negative openness needs
  expect_equal(max(read_band(rvt_resample(p, 8, method = "min"))), 0)
})

test_that("the clear-sky solar geometry matches published values", {
  # Kasten & Young air mass: 1 at the zenith, ~2 at 30 degrees, ~10 at 5
  expect_equal(.air_mass(90), 1, tolerance = 1e-3)
  expect_equal(.air_mass(30), 2, tolerance = 0.02)
  expect_equal(.air_mass(5), 10, tolerance = 0.1)
  # and it must grow monotonically as the sun drops
  a <- .air_mass(c(90, 60, 30, 10, 5))
  expect_true(all(diff(a) > 0))

  # Earth-Sun distance: perihelion in early January, about +/-3.3%
  e <- .extraterrestrial(as.Date(c("2025-01-03", "2025-07-04")))
  expect_gt(e[1], e[2])
  expect_equal(e[1] / e[2], 1.069, tolerance = 0.01)

  # diffuse rises with sun altitude and with haze
  expect_gt(.diffuse_horizontal(60, 1367, 3), .diffuse_horizontal(20, 1367, 3))
  expect_gt(.diffuse_horizontal(40, 1367, 6), .diffuse_horizontal(40, 1367, 2))
})

test_that("insolation on open flat ground is physically plausible", {
  flat <- make_raster(matrix(700, 40, 40))
  args <- list(flat, lat = 47.5, lon = 11.1, reach = 20, num_directions = 16)

  day <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), time_step = 10))))
  # a cloudless June day at 47 N delivers roughly 7-9 kWh/m2 on the horizontal
  expect_gt(mean(day), 7)
  expect_lt(mean(day), 10)

  # beam plus diffuse must equal global, cell for cell
  b <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), component = "beam"))))
  d <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), component = "diffuse"))))
  g <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), component = "global"))))
  expect_equal(b + d, g, tolerance = 1e-5)
  expect_gt(mean(b), mean(d))              # clear sky: beam dominates

  # "all" returns the same three as bands
  p <- do.call(rvt_insolation, c(args, list(dates = as.Date("2025-06-21"),
                                             component = "all")))
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  expect_equal(ds$getRasterCount(), 3L)
  ds$close()

  # haze cuts the beam and lifts the diffuse
  hazy_b <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), component = "beam", linke = 6))))
  hazy_d <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), component = "diffuse", linke = 6))))
  expect_lt(mean(hazy_b), mean(b))
  expect_gt(mean(hazy_d), mean(d))

  # units are just a rescaling of the same field
  mj <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), units = "mj"))))
  expect_equal(mj, g * 3.6, tolerance = 1e-4)
})

test_that("insolation responds to aspect, and more sharply in winter", {
  n <- 60
  fall <- tan(30 * pi / 180)
  south <- outer(seq(n, 1) * fall, rep(1, n)) + 700
  north <- outer(seq(1, n) * fall, rep(1, n)) + 700
  mid <- function(m) mean(m[25:35, 25:35])
  val <- function(z, d) mid(read_band(rvt_insolation(make_raster(z),
    dates = as.Date(d), lat = 47.5, lon = 11.1, reach = 20,
    num_directions = 16)))

  s6 <- val(south, "2025-06-21"); n6 <- val(north, "2025-06-21")
  s12 <- val(south, "2025-12-21"); n12 <- val(north, "2025-12-21")

  expect_gt(s6, n6)                        # south wins in both seasons
  expect_gt(s12, n12)
  expect_gt(s12 / n12, s6 / n6)            # but far more decisively in winter

  # in midwinter the north slope gets no direct sun at all, yet is not black:
  # the sky still lights it, and diffuse is not gated by the shadow test
  nb <- mid(read_band(rvt_insolation(make_raster(north),
    dates = as.Date("2025-12-21"), component = "beam", lat = 47.5, lon = 11.1,
    reach = 20, num_directions = 16)))
  expect_lt(nb, 0.05)
  expect_gt(n12, 0.2)
})

test_that("pyramid_method chooses how the coarse levels are built", {
  # max keeps a narrow obstruction, q3 discards it - that is the whole trade,
  # and which is right depends on whether the surface is bare earth or canopy
  z <- matrix(0, 64, 64)
  z[30:31, ] <- 50
  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = 64,
                           ysize = 64, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(557000, 1, 0, 5700000, 0, -1))
  ds$setProjection(gdalraster::srs_to_wkt("EPSG:25832"))
  ds$write(1L, 0, 0, 64, 64, as.vector(t(z)))
  ds$close()
  on.exit(unlink(p), add = TRUE)
  expect_equal(max(read_band(rvt_resample(p, 8, method = "max"))), 50)
  expect_equal(max(read_band(rvt_resample(p, 8, method = "q3"))), 0)

  # the argument reaches the levels: a rough surface coarsened by q3 must see
  # less obstruction than the same surface coarsened by max. Needs a CRS,
  # since these statistics are done with gdalwarp.
  set.seed(3)
  rough <- matrix(pmax(0, rnorm(400 * 400, 0, 6)), 400, 400)
  r <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = r, xsize = 400,
                           ysize = 400, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, r, read_only = FALSE)
  ds$setGeoTransform(c(557000, 1, 0, 5700000, 0, -1))
  ds$setProjection(gdalraster::srs_to_wkt("EPSG:25832"))
  ds$write(1L, 0, 0, 400, 400, as.vector(t(rough)))
  ds$close()
  on.exit(unlink(r), add = TRUE)

  mx <- read_band(rvt_svf(r, reach = 200, pyramid_method = "max"))
  q3 <- read_band(rvt_svf(r, reach = 200, pyramid_method = "q3"))
  expect_gt(mean(q3), mean(mx))

  # negative openness inverts the terrain, so its levels take the mirror
  # statistic; a method with no mirror is refused rather than silently wrong
  expect_equal(.mirror_method("max"), "min")
  expect_equal(.mirror_method("q3"), "q1")
  expect_equal(.mirror_method("med"), "med")
  expect_equal(.mirror_method("cubic"), "cubic")
  expect_error(.mirror_method("rms"), "does not commute")
  expect_error(rvt_openness_negative(dem, reach = 400, pyramid_method = "rms"),
               "does not commute")
})

test_that("rvt_fill() closes holes without touching measured ground", {
  z <- read_band(dem)
  f <- read_band(rvt_fill(dem))
  expect_equal(sum(is.na(z)), 393L)
  expect_equal(sum(is.na(f)), 0L)
  # cells that had data are returned exactly as they were
  expect_equal(f[!is.na(z)], z[!is.na(z)], tolerance = 1e-9)
  # and nothing is invented outside the range of the surrounding ground
  expect_gte(min(f), min(z, na.rm = TRUE))
  expect_lte(max(f), max(z, na.rm = TRUE))

  expect_error(rvt_fill(dem, max_distance = 0.5), "smaller than one cell")
})

test_that("smoothing removes noise but keeps a scarp sharp", {
  set.seed(5)
  n <- 200
  clean <- matrix(0, n, n)
  clean[, 101:n] <- 5                       # a one-cell, 5 m step
  noisy <- clean + matrix(rnorm(n * n, 0, 0.15), n, n)

  p <- tempfile(fileext = ".tif")
  ds <- gdalraster::create(format = "GTiff", dst_filename = p, xsize = n,
                           ysize = n, nbands = 1, dataType = "Float32")
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = FALSE)
  ds$setGeoTransform(c(557000, 1, 0, 5700000, 0, -1))
  ds$setProjection(gdalraster::srs_to_wkt("EPSG:25832"))
  ds$write(1L, 0, 0, n, n, as.vector(t(noisy)))
  ds$close()
  on.exit(unlink(p), add = TRUE)

  flat <- c(20:80, 120:180)
  noise_of <- function(m) stats::sd((m - clean)[, flat])
  # cells taken to cross the middle 80% of the step
  width_of <- function(m) {
    prof <- colMeans(m)
    lo <- min(prof); hi <- max(prof)
    sum(prof > lo + 0.1 * (hi - lo) & prof < lo + 0.9 * (hi - lo))
  }

  s <- read_band(rvt_smooth(p, radius = 5, norm_diff = 15, iterations = 3,
                             max_diff = 2))
  expect_lt(noise_of(s), noise_of(noisy) / 3)   # noise substantially down
  expect_equal(width_of(s), 0L)                 # step still one cell wide
  expect_equal(mean(s[, 130]) - mean(s[, 70]), 5, tolerance = 0.01)

  # this is the whole point: a Gaussian of comparable scale removes noise too,
  # but spreads that step over several cells
  g <- .gauss_kernel(2)
  gm <- noisy
  for (i in seq_len(n)) gm[i, ] <- stats::filter(gm[i, ], g, sides = 2)
  for (j in seq_len(n)) gm[, j] <- stats::filter(gm[, j], g, sides = 2)
  # stats::filter leaves NA at both margins, so measure on the interior
  prof <- colMeans(gm, na.rm = TRUE)
  prof <- prof[is.finite(prof)]
  lo <- min(prof); hi <- max(prof)
  expect_gt(sum(prof > lo + 0.1 * (hi - lo) & prof < lo + 0.9 * (hi - lo)), 3)

  # max_diff really is a cap
  s2 <- read_band(rvt_smooth(p, radius = 5, max_diff = 0.05))
  expect_lte(max(abs(s2 - noisy)), 0.05 + 1e-6)

  # and the result does not depend on how the raster was tiled
  expect_identical(read_band(rvt_smooth(p, radius = 5, tile_size = 64)),
                   read_band(rvt_smooth(p, radius = 5)))
})

test_that("smoothing leaves NoData alone", {
  z <- read_band(dem)
  s <- read_band(rvt_smooth(dem, radius = 3))
  expect_identical(is.na(s), is.na(z))
  expect_lte(max(abs(s - z), na.rm = TRUE), 0.5 + 1e-6)
  expect_error(rvt_smooth(dem, norm_diff = 0), "between 0 and 90")
  expect_error(rvt_smooth(dem, max_diff = -1), "positive")
})
