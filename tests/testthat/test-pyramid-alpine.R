# The multi-resolution search on real mountains. dtm1.tif has 71 m of relief
# over a kilometre, so distant terrain barely matters there and every
# reach-related claim would be validated where it cannot fail; the
# Partnachklamm tiles have 494 m over the same footprint, and give a genuine
# four-file mosaic rather than one raster cut into quadrants.
#
# Gitignored, so all three blocks skip when the tiles are not present.

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

