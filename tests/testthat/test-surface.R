# Surface shape: aspect, TPI/DEV, geomorphons, curvature and
# Laplacian-of-Gaussian. Split out of test-reference.R, which had grown to
# half the suite in one file.
#
# The curvature blocks check against closed-form values on analytic surfaces
# rather than against another implementation, because the naming across
# ArcGIS, SAGA, GRASS and WhiteboxTools is a swamp - see CLAUDE.md. Shared
# helpers (read_band, make_raster, dem) live in helper-rvtr.R.

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

test_that("a windowed curvature fit is exact on quadratic surfaces", {
  # The fit recovers any quadratic exactly, at *any* radius - that is what
  # makes these surfaces the right check: unlike the 3x3 kernel, where a wider
  # window would be an approximation, here there is nothing to approximate, so
  # the tolerance can be tight rather than the 1e-2 the block above needs.
  # Values are read at the centre cell alone, because the analytic curvature of
  # a tilted paraboloid varies across the raster and a median over any patch
  # would not equal it.
  n <- 121L; ctr <- 61L
  x <- col(matrix(0, n, n)) - ctr
  y <- row(matrix(0, n, n)) - ctr
  at_ctr <- function(p, ty, r) read_band(rvt_curvature(p, type = ty,
                                                        radius = r))[ctr, ctr]

  a <- 0.01; tilt <- 0.3; q <- 1 + tilt^2
  pb <- make_raster(matrix(100 + a * (x^2 + y^2) + tilt * x, n, n))
  H <- -((1 + tilt^2) * (2 * a) + (2 * a)) / (2 * q^1.5)
  K <- (2 * a * 2 * a) / q^2
  for (r in c(2, 5, 10)) {
    expect_equal(at_ctr(pb, "mean", r), H, tolerance = 1e-4, label = paste("mean r", r))
    expect_equal(at_ctr(pb, "gaussian", r), K, tolerance = 1e-4,
                 label = paste("gaussian r", r))
  }

  # the saddle is the only one of these that exercises the mixed partial
  cc <- 0.01
  sd_ <- make_raster(matrix(100 + cc * x * y, n, n))
  expect_lt(abs(at_ctr(sd_, "mean", 5)), 1e-9)
  expect_equal(at_ctr(sd_, "gaussian", 5), -cc^2, tolerance = 1e-4)

  # a plane has no curvature however wide the window
  pl <- make_raster(matrix(100 + 0.2 * x + 0.1 * y, n, n))
  for (ty in c("profile", "plan", "mean", "total", "gaussian"))
    expect_lt(abs(at_ctr(pl, ty, 5)), 1e-9, label = ty)

  unlink(c(pb, sd_, pl))
})

test_that("Shary and Koenderink curvatures have their analytic values", {
  # These four are arithmetic on the principal curvatures, so the shapes below
  # pin them exactly: a paraboloid is umbilic at its centre (both principals
  # equal) and a saddle is the opposite extreme (equal and opposite).
  n <- 121L; ctr <- 61L
  x <- col(matrix(0, n, n)) - ctr
  y <- row(matrix(0, n, n)) - ctr
  at_ctr <- function(p, ty) read_band(rvt_curvature(p, type = ty))[ctr, ctr]

  # Tolerance on the *magnitudes* is 1e-3, not tighter: the surfaces are stored
  # Float32, so elevations near 100 quantise at ~8e-6 and a second difference
  # over 1 m cells turns that into a few parts in 10,000 of a 0.02 curvature
  # (measured: 2.1e-4). The ratios below are unaffected, because the same
  # quantisation cancels top and bottom, so they keep a tight tolerance.

  # bowl z = a(x^2+y^2): umbilic, and concave, so convex-positive makes it a
  # pit. kmax = kmin = -2a.
  a <- 0.01
  bowl <- make_raster(matrix(100 + a * (x^2 + y^2), n, n))
  expect_equal(at_ctr(bowl, "unsphericity"), 0, tolerance = 1e-6)
  expect_equal(at_ctr(bowl, "casorati"), 2 * a, tolerance = 1e-3)
  expect_equal(at_ctr(bowl, "shape_index"), -1, tolerance = 1e-4)
  # and the dome is its mirror in every one of them
  dome <- make_raster(matrix(100 - a * (x^2 + y^2), n, n))
  expect_equal(at_ctr(dome, "shape_index"), 1, tolerance = 1e-4)
  expect_equal(at_ctr(dome, "casorati"), 2 * a, tolerance = 1e-3)

  # saddle z = c.x.y: kmax = -kmin = c, so a perfect saddle reads 0
  cc <- 0.01
  sd_ <- make_raster(matrix(100 + cc * x * y, n, n))
  expect_equal(at_ctr(sd_, "unsphericity"), cc, tolerance = 1e-3)
  expect_equal(at_ctr(sd_, "casorati"), cc, tolerance = 1e-3)
  expect_equal(at_ctr(sd_, "shape_index"), 0, tolerance = 1e-6)

  # a plane bends in no direction at all
  pl <- make_raster(matrix(100 + 0.2 * x + 0.1 * y, n, n))
  expect_equal(at_ctr(pl, "unsphericity"), 0, tolerance = 1e-9)
  expect_equal(at_ctr(pl, "casorati"), 0, tolerance = 1e-9)

  # unsphericity and casorati are magnitudes; on real terrain neither may go
  # negative, and the shape index must stay inside its bounds
  for (ty in c("unsphericity", "casorati")) {
    v <- read_band(fx_curvature(ty))
    expect_true(all(v[!is.na(v)] >= 0), label = ty)
  }
  si <- read_band(fx_curvature("shape_index"))
  expect_true(all(abs(si[!is.na(si)]) <= 1))

  unlink(c(bowl, dome, sd_, pl))
})

test_that("the landform signatures in ?rvt_curvature hold", {
  # The help describes each curvature by what it says about the terrain. These
  # are those sentences turned into assertions, so the prose cannot drift from
  # the kernel - it already had: it claimed a ridge crest has "a large negative
  # minimal and a near-zero maximal", which is a valley, exactly backwards.
  n <- 121L; ctr <- 61L
  x <- col(matrix(0, n, n)) - ctr
  y <- row(matrix(0, n, n)) - ctr
  sig <- function(z, ty) {
    p <- make_raster(matrix(z, n, n))
    on.exit(unlink(p), add = TRUE)
    read_band(rvt_curvature(p, type = ty))[ctr, ctr]
  }
  ridge  <- 100 - 0.02 * x^2
  valley <- 100 + 0.02 * x^2
  dome   <- 100 - 0.01 * (x^2 + y^2)
  pit    <- 100 + 0.01 * (x^2 + y^2)
  saddle <- 100 + 0.01 * x * y

  # "a ridge bends hard across its crest and barely at all along it"
  expect_gt(sig(ridge, "maximal"), 0.03)
  expect_equal(sig(ridge, "minimal"), 0, tolerance = 1e-6)
  # "a valley is the mirror of that"
  expect_lt(sig(valley, "minimal"), -0.03)
  expect_equal(sig(valley, "maximal"), 0, tolerance = 1e-6)
  # "a dome has both positive, a pit both negative"
  expect_gt(sig(dome, "minimal"), 0)
  expect_lt(sig(pit, "maximal"), 0)
  # gaussian: same sign every way is dome/bowl, opposite ways is a saddle
  expect_gt(sig(dome, "gaussian"), 0)
  expect_gt(sig(pit, "gaussian"), 0)
  expect_lt(sig(saddle, "gaussian"), 0)
  # unsphericity: zero where the ground curves alike in every direction
  expect_equal(sig(dome, "unsphericity"), 0, tolerance = 1e-6)
  expect_gt(sig(ridge, "unsphericity"), 0.01)
  expect_gt(sig(saddle, "unsphericity"), 0.005)
  # the documented shape_index scale: pit -1, valley -0.5, saddle 0,
  # ridge +0.5, peak +1
  expect_equal(sig(pit, "shape_index"),    -1,   tolerance = 1e-4)
  expect_equal(sig(valley, "shape_index"), -0.5, tolerance = 1e-4)
  expect_equal(sig(saddle, "shape_index"),  0,   tolerance = 1e-6)
  expect_equal(sig(ridge, "shape_index"),   0.5, tolerance = 1e-4)
  expect_equal(sig(dome, "shape_index"),    1,   tolerance = 1e-4)
})

test_that("twisting and rotor match their analytic identities", {
  # The normalisation is the whole difficulty with these two, so they are
  # pinned by identities rather than by remembered formulae.
  n <- 81L; ctr <- 41L
  x <- col(matrix(0, n, n)) - ctr
  y <- row(matrix(0, n, n)) - ctr          # +y is south, as the kernel reads it
  at_ctr <- function(z, ty) {
    p <- make_raster(matrix(z, n, n))
    on.exit(unlink(p), add = TRUE)
    read_band(rvt_curvature(p, type = ty))[ctr, ctr]
  }

  # A ramp that twists: z = -x - a.x.y. At the origin the gradient is (-1, 0),
  # so downhill is due east, and the slope steepens towards the south.
  a <- 0.1
  ramp <- -x - a * x * y
  # twisting *is* d(slope angle)/d(arc length) along the contour. Here the
  # contour is the line x = 0, on which z is identically 0, so that rate is
  # a/(1+1) = a/2 in closed form.
  expect_equal(at_ctr(ramp, "twisting"), a / 2, tolerance = 1e-2)
  # rotor is the same numerator with the projected normalisation, = a here
  expect_equal(at_ctr(ramp, "rotor"), a, tolerance = 1e-2)
  # and both are positive, meaning "to the right looking downhill" - the
  # descent direction really does swing southward, which is right of east
  expect_gt(at_ctr(ramp, "twisting"), 0)
  expect_gt(at_ctr(ramp, "rotor"), 0)
  # mirroring north for south must flip both and nothing else
  expect_equal(at_ctr(-x + a * x * y, "twisting"), -a / 2, tolerance = 1e-2)
  expect_equal(at_ctr(-x + a * x * y, "rotor"), -a, tolerance = 1e-2)

  # Surfaces whose contours are already lines of curvature cannot twist:
  # a cone/dome (circular contours) and a cylinder ridge (straight ones).
  expect_equal(at_ctr(100 - 0.01 * (x^2 + y^2), "twisting"), 0, tolerance = 1e-6)
  expect_equal(at_ctr(100 - 0.02 * x^2 + 0.3 * y, "twisting"), 0, tolerance = 1e-6)
  # a plane neither twists nor bends its flow lines
  expect_equal(at_ctr(100 + 0.2 * x + 0.1 * y, "twisting"), 0, tolerance = 1e-9)
  expect_equal(at_ctr(100 + 0.2 * x + 0.1 * y, "rotor"), 0, tolerance = 1e-9)
})

test_that("twisting never exceeds unsphericity", {
  # In the principal frame the geodesic torsion is ((k2-k1)/2).sin(2.theta),
  # so |twisting| <= unsphericity everywhere, with equality where the contour
  # bisects the principal directions. This is the check that would catch a
  # wrong denominator on real terrain, which the synthetic shapes above cannot.
  tw <- read_band(fx_curvature("twisting"))
  un <- read_band(fx_curvature("unsphericity"))
  ok <- !is.na(tw) & !is.na(un)
  # a hair of slack for the two rasters being stored Float32 independently
  expect_true(all(abs(tw[ok]) <= un[ok] + 1e-6))
  # the bound is tight, not vacuous: somewhere it is very nearly attained
  expect_gt(max((abs(tw[ok]) / un[ok])[un[ok] > 1e-3]), 0.95)
})

test_that("difference curvature is the split mean curvature hides", {
  # kd = (profile - tangential)/2 by definition, so the check is that the
  # kernel's own profile and tangential reproduce it exactly - a consistency
  # test, which is the right kind here because that *is* the definition.
  g <- function(ty, radius = NULL) {
    m <- read_band(fx_curvature(ty, radius))[100:900, 100:900]
    m[!is.na(m)]
  }
  expect_equal(g("difference"), (g("profile") - g("tangential")) / 2,
               tolerance = 1e-6)
  # and it must hold on the windowed kernel too, not just the 3x3 one
  expect_equal(g("difference", radius = 5),
               (g("profile", radius = 5) - g("tangential", radius = 5)) / 2,
               tolerance = 1e-6)
})

test_that("the curvature default is a different estimator from the fit", {
  # Omitting `radius` runs the 3x3 second difference (what ArcGIS and
  # WhiteboxTools return); giving one fits a quadratic, which at r = 1 averages
  # the bending over all three rows instead of taking the middle row alone.
  # An isolated one-cell bump separates them exactly: the fit sees two flat
  # rows beside the spike and reports a third of the second difference.
  z <- matrix(0, 41L, 41L); z[21L, 21L] <- 1
  bump <- make_raster(z)
  on.exit(unlink(bump), add = TRUE)
  # total = zxx^2 + 2 zxy^2 + zyy^2, and zxy is 0 here, so sqrt(total/2) is |zxx|
  zxx <- function(...) sqrt(read_band(rvt_curvature(bump, type = "total",
                                                     ...))[21L, 21L] / 2)
  expect_equal(zxx(), 2)                      # left + right - 2*centre
  expect_equal(zxx(radius = 1), 2 / 3)        # averaged over three rows

  # on real terrain the two must therefore disagree, and a wider window must
  # damp the fine detail the finest scale is dominated by
  prof <- function(radius = NULL) {
    m <- read_band(fx_curvature("profile", radius))[100:900, 100:900]
    m[!is.na(m)]
  }
  p0 <- prof(); p1 <- prof(radius = 1); p5 <- prof(radius = 5)
  expect_false(isTRUE(all.equal(p0, p1)))
  expect_gt(sd(p0), sd(p1))
  expect_gt(sd(p1), sd(p5))
})

test_that("windowed curvature tiles invariantly and keeps the NoData mask", {
  # A new overlap path - the default kernel pads by one cell, the fit by the
  # whole radius - so the package's core guarantee needs re-checking here.
  whole <- read_band(fx_curvature("profile", 5))
  tiled <- read_band(rvt_curvature(dem, radius = 5, tile_size = 137))
  expect_identical(whole, tiled)

  # NoData inside the window contributes the centre cell's own elevation, so a
  # hole must not erode a ring of cells around itself the way it once did for
  # slope (see deviation 4 in CLAUDE.md)
  expect_identical(is.na(whole), is.na(read_band(dem)))
})

test_that("Laplacian-of-Gaussian is flat on a plane and signed across an edge", {
  n <- 81L
  x <- col(matrix(0, n, n)) - 41
  pl <- make_raster(matrix(100 + 0.2 * x, n, n))
  expect_lt(max(abs(read_band(rvt_log(pl, sigma = 2))[30:52, 30:52])), 1e-6)

  step <- matrix(100, n, n); step[, 42:n] <- 101
  sp <- make_raster(step)
  lg <- read_band(rvt_log(sp, sigma = 2))[41, ]
  # a step produces the classic antisymmetric pair, summing to ~0. Convex is
  # positive, as for rvt_curvature(): the lip on the high side is positive,
  # the foot on the low side negative - the negative of the raw Laplacian
  expect_lt(min(lg[35:41]), 0)
  expect_gt(max(lg[42:48]), 0)
  expect_equal(sum(lg[35:48]), 0, tolerance = 1e-6)

  # and the same sign as curvature on the same edge
  cv <- read_band(rvt_curvature(sp))[41, ]
  expect_equal(sign(lg[c(39, 44)]), sign(cv[c(41, 42)]))

  unlink(c(pl, sp))
})

