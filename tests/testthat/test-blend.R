# Layer blending. Two things carry most of the weight here:
#
#  - every blend mode against its closed form, because thirteen near-identical
#    algebraic formulas are exactly where a transcription slip hides and no
#    amount of looking at output would reveal it;
#  - the engine reproducing .vat_stack(), which is already validated against
#    rvt-py. That pins the whole pipeline - normalisation, mode, opacity,
#    ordering - to a known-good result without touching rvt_vat() itself.

dem <- system.file("extdata", "dtm1.tif", package = "rvtr")

read_b <- function(p, b = 1L) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  m <- matrix(as.numeric(ds$read(b, 0, 0, nx, ny, nx, ny)),
              nrow = ny, ncol = nx, byrow = TRUE)
  nd <- ds$getNoDataValue(b)
  if (!is.na(nd)) m[m == nd] <- NA
  m
}
n_bands <- function(p) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  ds$getRasterCount()
}

test_that("every blend mode matches its closed form", {
  # endpoints and midpoints, including the 0 and 1 that make dodge and burn
  # divide by zero if they are not guarded
  a <- c(0, 0.25, 0.5, 0.75, 1, 0.3, 0.8)
  b <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.9, 0.1)
  m <- rvtr:::.blend_modes

  expect_equal(m$normal(a, b),     a)
  expect_equal(m$lighten(a, b),    pmax(a, b))
  expect_equal(m$darken(a, b),     pmin(a, b))
  expect_equal(m$multiply(a, b),   a * b)
  expect_equal(m$screen(a, b),     1 - (1 - a) * (1 - b))
  expect_equal(m$addition(a, b),   pmin(1, a + b))
  expect_equal(m$subtract(a, b),   pmax(0, b - a))
  expect_equal(m$difference(a, b), abs(a - b))
  expect_equal(m$overlay(a, b),    ifelse(b <= 0.5, 2 * a * b,
                                           1 - 2 * (1 - a) * (1 - b)))
  expect_equal(m$hard_light(a, b), ifelse(a <= 0.5, 2 * a * b,
                                           1 - 2 * (1 - a) * (1 - b)))
  expect_equal(m$dodge(a, b),      pmin(1, ifelse(a >= 1, 1, b / (1 - a))))
  expect_equal(m$burn(a, b),       pmax(0, ifelse(a <= 0, 0, 1 - (1 - b) / a)))
  expect_equal(m$soft_light(a, b), {
    d <- ifelse(b <= 0.25, ((16 * b - 12) * b + 4) * b, sqrt(b))
    ifelse(a <= 0.5, b - (1 - 2 * a) * b * (1 - b), b + (2 * a - 1) * (d - b))
  })

  # structural identities, which catch a swapped argument that the formulas
  # above would not: screen is multiply on the inverted inputs, and hard light
  # is overlay with the two layers exchanged
  expect_equal(m$screen(a, b), 1 - m$multiply(1 - a, 1 - b))
  expect_equal(m$hard_light(a, b), m$overlay(b, a))

  # nothing may leave 0-1, or it would clip unpredictably once written
  for (nm in rvt_blend_modes) {
    v <- m[[nm]](a, b)
    expect_true(all(v >= 0 & v <= 1), label = nm)
  }
})

test_that("the blend engine reproduces the VAT stack", {
  # This is the load-bearing test. .vat_stack() is checked against rvt-py, so
  # matching it exactly means the engine's normalisation, mode arithmetic,
  # opacity and layer ordering are all right at once.
  hs   <- rvt_hillshade(dem, sun_elevation = 35)
  slp  <- rvt_slope(dem)
  svf  <- rvt_svf(dem, reach = 10)
  opns <- rvt_openness(dem, reach = 10)
  on.exit(unlink(c(hs, slp, svf, opns)), add = TRUE)

  p <- rvt_preset_general
  got <- rvt_stack(hs, range = c(0, 1)) |>
    rvt_blend(slp,  "normal",   0.5,  range = p$slope, invert = TRUE) |>
    rvt_blend(opns, "overlay",  0.5,  range = p$opns) |>
    rvt_blend(svf,  "multiply", 0.25, range = p$svf) |>
    rvt_render()
  want <- rvtr:::.vat_stack(read_b(hs), read_b(slp), read_b(opns),
                            read_b(svf), p)

  a <- read_b(got)
  ok <- !is.na(a) & !is.na(want)
  # 1e-6 is the Float32 output floor, not slack in the arithmetic
  expect_lt(max(abs(a[ok] - want[ok])), 1e-6)
  expect_identical(is.na(a), is.na(want))
})

test_that("blending is tile-invariant and preserves NoData", {
  hs  <- rvt_hillshade(dem)
  svf <- rvt_svf(dem, reach = 10)
  on.exit(unlink(c(hs, svf)), add = TRUE)

  st <- hs |> rvt_blend(svf, "multiply", 0.25, range = c(0.7, 1))
  # ranges are frozen before the first tile, which is the whole reason this
  # holds - a per-tile stretch would give the same ground different values
  expect_identical(read_b(rvt_render(st)),
                   read_b(rvt_render(st, tile_size = 137)))
  # a blend of a hole is a hole, and only a hole
  expect_identical(is.na(read_b(rvt_render(st))), is.na(read_b(dem)))
})

test_that("rvt_blend() builds a recipe and writes nothing until rendered", {
  hs  <- rvt_hillshade(dem)
  svf <- rvt_svf(dem, reach = 10)
  on.exit(unlink(c(hs, svf)), add = TRUE)

  d <- fs::path_temp()
  before <- length(fs::dir_ls(d, glob = "*.tif"))
  st <- hs |> rvt_blend(svf, "multiply", 0.25) |> rvt_blend(svf, "screen", 0.1)
  expect_s3_class(st, "rvt_stack")
  expect_length(st$layers, 3L)
  # deferring the write is the entire justification for returning an object
  expect_equal(length(fs::dir_ls(d, glob = "*.tif")), before)

  rvt_render(st)
  expect_equal(length(fs::dir_ls(d, glob = "*.tif")), before + 1L)

  # a bare path and rvt_stack() must be the same starting point
  expect_identical(read_b(rvt_render(rvt_blend(hs, svf, "multiply", 0.25))),
                   read_b(rvt_render(rvt_blend(rvt_stack(hs), svf,
                                                "multiply", 0.25))))
})

test_that("opacity 0 leaves the background untouched, for every mode", {
  hs  <- rvt_hillshade(dem)
  svf <- rvt_svf(dem, reach = 10)
  on.exit(unlink(c(hs, svf)), add = TRUE)
  base <- read_b(rvt_render(rvt_stack(hs)))
  for (nm in rvt_blend_modes)
    expect_equal(read_b(rvt_render(rvt_blend(hs, svf, nm, opacity = 0))),
                 base, tolerance = 1e-6, label = nm)
})

test_that("a three-band layer broadcasts against a single-band one", {
  hs  <- rvt_hillshade(dem)
  rgb <- rvt_mstp(dem)
  on.exit(unlink(c(hs, rgb)), add = TRUE)
  expect_equal(n_bands(rgb), 3L)

  out <- rgb |> rvt_blend(hs, "multiply", opacity = 1) |> rvt_render()
  expect_equal(n_bands(out), 3L)

  # one range spans all three bands, so the colour balance is preserved; a
  # per-band stretch would silently change the hue
  rs <- vapply(1:3, function(b) rvt_range(rgb, band = b), numeric(2))
  rng <- c(min(rs[1, ]), max(rs[2, ]))
  a <- rvtr:::.norm_lin(read_b(hs), rvt_range(hs)[1], rvt_range(hs)[2])
  for (k in 1:3) {
    want <- a * rvtr:::.norm_lin(read_b(rgb, k), rng[1], rng[2])
    got <- read_b(out, k)
    ok <- !is.na(want) & !is.na(got)
    expect_lt(max(abs(want[ok] - got[ok])), 1e-6)
  }

  # and the other way round: a 3-band layer over a 1-band background
  expect_equal(n_bands(rvt_render(rvt_blend(hs, rgb, "screen", 0.5))), 3L)
})

test_that("misaligned layers are refused, not silently blended", {
  hs <- rvt_hillshade(dem)
  small <- tempfile(fileext = ".tif")
  gdalraster::translate(dem, small, quiet = TRUE,
                         cl_arg = c("-srcwin", "0", "0", "500", "500"))
  # same size and cell size as the hillshade, but shifted 10 m east: the old
  # check compared only size and resolution, and let this through
  shifted <- tempfile(fileext = ".tif")
  gdalraster::translate(hs, shifted, quiet = TRUE)
  ds <- methods::new(gdalraster::GDALRaster, shifted, read_only = FALSE)
  gt <- ds$getGeoTransform(); gt[1] <- gt[1] + 10; ds$setGeoTransform(gt); ds$close()
  # 2.5 m against 1 m: not a whole factor
  r25 <- rvt_resample(dem, 2.5)
  on.exit(unlink(c(hs, small, shifted, r25)), add = TRUE)

  expect_error(rvt_blend(hs, small, "multiply"), "same extent")
  expect_error(rvt_blend(hs, shifted, "multiply"), "same extent")
  expect_error(rvt_blend(hs, r25, "multiply"), "whole")
})

test_that("layers at whole-factor resolutions render on the finest grid", {
  hs <- rvt_hillshade(dem)
  svf4 <- rvt_svf(rvt_resample(dem, 4), reach = 16)       # 250 x 250 at 4 m
  on.exit(unlink(c(hs, svf4)), add = TRUE)

  # whichever way round, the output is the fine grid
  a <- rvt_render(rvt_blend(hs, svf4, "multiply", opacity = 0.5))
  b <- rvt_render(rvt_blend(svf4, hs, "multiply", opacity = 0.5))
  on.exit(unlink(c(a, b)), add = TRUE)
  expect_equal(dim(read_b(a)), c(1000L, 1000L))
  expect_equal(dim(read_b(b)), c(1000L, 1000L))
  expect_output(print(rvt_blend(svf4, hs)), "1000 x 1000, 1 m")

  # the coarse layer arrives smoothly, not as 4 x 4 blocks: bilinear, so
  # neighbouring fine cells inside one coarse cell differ
  up <- read_b(rvt_render(rvt_stack(svf4) |> rvt_blend(hs, opacity = 0)))
  expect_gt(mean(up[2:999, 2] != up[2:999, 3]), 0.5)

  # at opacity 0 the coarse layer changes nothing where it has data (NoData in
  # any layer is NoData out, whatever its opacity)
  fine_only <- read_b(rvt_render(rvt_stack(hs)))
  zero <- read_b(rvt_render(rvt_blend(hs, svf4, "multiply", opacity = 0)))
  ok <- !is.na(zero)
  expect_gt(mean(ok), 0.95)
  expect_equal(zero[ok], fine_only[ok])

  # and the render still cannot depend on tile_size
  whole  <- read_b(a)
  pieces <- read_b(rvt_render(rvt_blend(hs, svf4, "multiply", opacity = 0.5),
                              tile_size = 137))
  expect_identical(whole, pieces)
})

test_that("rvt_range reads a frozen stretch off the raster", {
  svf <- rvt_svf(dem, reach = 10)
  on.exit(unlink(svf), add = TRUE)
  full <- rvt_range(svf)
  cut <- rvt_range(svf, pct = c(2, 98))
  expect_length(full, 2L)
  expect_lt(full[1], full[2])
  # a percentile cut sits strictly inside the full range, which is the point
  expect_gt(cut[1], full[1])
  expect_lt(cut[2], full[2])
  expect_error(rvt_range(svf, pct = c(98, 2)), "increasing")
})
