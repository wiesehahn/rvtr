# Preprocessing: rvt_fill() (GDAL inverse-distance interpolation) and
# rvt_smooth() (Sun, Rosin, Martin and Langbein 2007, the algorithm behind
# WhiteboxTools FeaturePreservingSmoothing).
#
# The trap rvt_smooth() had to be rewritten for: the elevation update must
# use the immediate 3x3 neighbours only. Letting it reach the full `radius`
# window smears the very edge the normal filter went to trouble to keep.

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
