# rvt_relight(). The output is for viewing, so these check behaviour rather
# than exact values. Network tests run on a small crop of one LGLN tile beside
# Burg Hardenberg to keep the suite fast.

relight_pt <- c(9.9464, 51.6317)
relight_rgb <- function() rvt_data_lgln(relight_pt, "rgb", year = 2025, res = 1)
relight_dsm <- function() rvt_data_lgln(relight_pt, "dsm")

relight_crop <- function(path, size = 300L) {
  out <- tempfile(fileext = ".tif")
  gdalraster::translate(path, out, quiet = TRUE,
                        cl_arg = c("-srcwin", "350", "350", size, size))
  out
}
relight_read <- function(p) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  list(nx = nx, ny = ny, nb = ds$getRasterCount(), gt = ds$getGeoTransform(),
       structure = ds$getMetadata(0, "IMAGE_STRUCTURE"),
       bands = lapply(seq_len(ds$getRasterCount()), function(b)
         as.numeric(ds$read(b, 0, 0, nx, ny, nx, ny))))
}

test_that("styles resolve, and explicit settings override them", {
  s <- rvtr:::.relight_params("default", list())
  expect_true(all(rvtr:::.relight_args %in% names(s)))

  w <- rvtr:::.relight_params("winter", list())
  expect_equal(w$sun_elevation, 18)
  expect_equal(w$bounce, rvt_relight_styles$default$bounce)   # inherited

  o <- rvtr:::.relight_params("winter", list(sun_elevation = 10, contrast = NULL))
  expect_equal(o$sun_elevation, 10)                           # override wins
  expect_equal(o$contrast, rvt_relight_styles$default$contrast)

  # every style must resolve to a usable set of settings
  for (nm in names(rvt_relight_styles))
    expect_silent(rvtr:::.relight_params(nm, list()))
})

test_that("bad styles and settings fail clearly, before any work", {
  expect_error(rvt_relight("a.tif", "b.tif", style = "sunset"),
               "Unknown style. Choose one of: default")
  expect_error(rvt_relight("a.tif", "b.tif", sun_elevation = 1, softness = 4),
               "between 0 and 90")
  expect_error(rvt_relight("a.tif", "b.tif", quality = 0), "quality")
})

test_that("every style renders a WebP COG on the orthophoto's grid", {
  skip_if_offline()
  rgb <- relight_crop(relight_rgb())
  dsm <- relight_crop(relight_dsm())
  on.exit(unlink(c(rgb, dsm)), add = TRUE)
  ref <- relight_read(rgb)

  for (nm in names(rvt_relight_styles)) {
    r <- relight_read(rvt_relight(rgb, dsm, style = nm))
    expect_equal(c(r$nx, r$ny, r$nb), c(ref$nx, ref$ny, 3), label = nm)
    expect_equal(r$gt, ref$gt, label = nm)
    expect_true(all(c("LAYOUT=COG", "COMPRESSION=WEBP") %in% r$structure),
                label = nm)
  }
})

test_that("neutral light gives the orthophoto back", {
  skip_if_offline()
  rgb <- relight_crop(relight_rgb())
  dsm <- relight_crop(relight_dsm())
  on.exit(unlink(c(rgb, dsm)), add = TRUE)

  # no sun, white sky, no occlusion: the light is exactly 1 everywhere
  out <- rvt_relight(rgb, dsm, saturation = 1, contrast = 1,
                     sun_strength = 0, sky_strength = 1, sky_color = "#FFFFFF",
                     occlusion = 0, quality = 100)
  a <- relight_read(rgb)$bands; b <- relight_read(out)$bands
  # WebP is lossy even at quality 100, so allow a small average error
  for (k in 1:3) expect_lt(mean(abs(a[[k]] - b[[k]])), 3)
})

test_that("tile borders leave no seams", {
  skip_if_offline()
  rgb <- relight_crop(relight_rgb())
  dsm <- relight_crop(relight_dsm())
  on.exit(unlink(c(rgb, dsm)), add = TRUE)

  # the same image rendered whole and in small 137-pixel pieces, an odd size
  # so the piece borders land in awkward places
  whole  <- relight_read(rvt_relight(rgb, dsm))$bands
  pieces <- relight_read(rvt_relight(rgb, dsm, tile_size = 137))$bands
  for (k in 1:3) expect_lte(max(abs(whole[[k]] - pieces[[k]])), 1)
})

test_that("an orthophoto on a different grid is refused", {
  skip_if_offline()
  rgb <- relight_crop(relight_rgb(), 300L)
  dsm <- relight_crop(relight_dsm(), 200L)
  on.exit(unlink(c(rgb, dsm)), add = TRUE)
  expect_error(rvt_relight(rgb, dsm), "rvt_resample")
})
