# Display images: WebP and JPEG written straight from the tiled pass, by
# rvt_render(), rvt_vat(), rvt_relight() and rvt_image().

dem <- system.file("extdata", "dtm1.tif", package = "rvtr")

img_read <- function(p) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  list(nx = nx, ny = ny, driver = ds$getDriverShortName(),
       bands = lapply(seq_len(ds$getRasterCount()), function(b)
         matrix(as.numeric(ds$read(b, 0, 0, nx, ny, nx, ny)), ny, nx, byrow = TRUE)))
}
cog_read <- function(p) {
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  m <- matrix(as.numeric(ds$read(1, 0, 0, nx, ny, nx, ny)), ny, nx, byrow = TRUE)
  nd <- ds$getNoDataValue(1)
  if (!is.na(nd)) m[m == nd] <- NA
  m
}
img_dir <- function() { d <- tempfile("img"); dir.create(d); d }

test_that("a blend renders straight to WebP and JPEG, with no sidecar files", {
  d <- img_dir(); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  st <- rvt_blend(rvt_hillshade(dem), rvt_svf(dem), "multiply", opacity = 0.5)
  cog <- cog_read(rvt_render(st))

  w <- img_read(rvt_render(st, file.path(d, "a.webp")))
  j <- img_read(rvt_render(st, file.path(d, "a.jpg")))
  expect_equal(c(w$nx, w$ny, length(w$bands)), c(1000, 1000, 4))
  expect_equal(w$driver, "WEBP")
  expect_equal(c(j$nx, j$ny, length(j$bands)), c(1000, 1000, 1))
  expect_equal(j$driver, "JPEG")
  # plain pictures: nothing but the two images in the folder
  expect_setequal(list.files(d), c("a.webp", "a.jpg"))

  # alpha is exact - lossless in WebP - and marks NoData precisely
  expect_identical(w$bands[[4]] == 0, is.na(cog))
  # lossy, but close to the 0-1 render scaled to 0-255
  ok <- !is.na(cog)
  expect_lt(mean(abs(w$bands[[1]][ok] - cog[ok] * 255)), 3)
  expect_lt(mean(abs(j$bands[[1]][ok] - cog[ok] * 255)), 4)
})

test_that("tiling changes nothing in the picture", {
  d <- img_dir(); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  st <- rvt_blend(rvt_hillshade(dem), rvt_svf(dem), "multiply", opacity = 0.5)
  # the encoder sees identical pixels either way, so even lossy output matches
  w <- img_read(rvt_render(st, file.path(d, "a.webp")))
  pieces <- img_read(rvt_render(st, file.path(d, "b.webp"), tile_size = 137))
  expect_identical(pieces$bands, w$bands)
})

test_that("metric functions refuse an image path and point to rvt_image()", {
  expect_error(rvt_svf(dem, tempfile(fileext = ".webp")), "rvt_image")
  expect_error(rvt_slope(dem, tempfile(fileext = ".jpg")), "rvt_image")
  expect_error(rvt_image(dem, tempfile(fileext = ".tif")), "webp")
  expect_error(rvt_render(rvt_stack(dem), tempfile(fileext = ".webp"), quality = 0),
               "quality")
})

test_that("the WebP size limit stops the call before any tile is computed", {
  big <- tempfile(fileext = ".vrt")
  gdalraster::translate(dem, big, quiet = TRUE,
                        cl_arg = c("-of", "VRT", "-outsize", "17000", "10"))
  on.exit(unlink(big), add = TRUE)
  t0 <- Sys.time()
  expect_error(rvt_image(big, tempfile(fileext = ".webp")), "16,383")
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 5)
  # JPEG allows it
  expect_silent(rvtr:::.image_check_size(17000, 10, "JPEG"))
})

test_that("rvt_image stretches over the range and colours through the palette", {
  d <- img_dir(); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  svf <- rvt_svf(dem)
  v <- cog_read(svf)

  # Grays over c(0.7, 1): the picture is lossy, so close rather than exact
  g <- img_read(rvt_image(svf, file.path(d, "g.webp"), range = c(0.7, 1)))
  expect_equal(length(g$bands), 4L)
  tab <- rvtr:::.col_table("Grays")
  ok <- !is.na(v)
  s <- pmin(pmax((v[ok] - 0.7) / 0.3, 0), 1)
  want <- round(tab[1, pmin(floor(s * 256), 255) + 1] * 255)
  expect_lt(mean(abs(g$bands[[1]][ok] - want)), 3)
  # the stretch ends land on the palette's end colours
  expect_lt(abs(mean(g$bands[[1]][ok][v[ok] >= 1]) - tab[1, 256] * 255), 3)

  # a named palette really is used: red, green and blue differ
  vir <- img_read(rvt_image(svf, file.path(d, "v.webp"), "Viridis"))
  expect_gt(mean(abs(vir$bands[[1]] - vir$bands[[2]])), 20)

  # a three-band raster is drawn as RGB, not through the palette: grey stays grey
  rgb <- tempfile(fileext = ".tif")
  gdalraster::translate(svf, rgb, quiet = TRUE, cl_arg = c("-b", "1", "-b", "1", "-b", "1"))
  on.exit(unlink(rgb), add = TRUE)
  r <- img_read(rvt_image(rgb, file.path(d, "rgb.webp"), "Viridis"))
  expect_lt(mean(abs(r$bands[[1]] - r$bands[[2]])), 1)
})

test_that("rvt_vat and rvt_image of a stack write pictures too", {
  d <- img_dir(); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  v <- img_read(rvt_vat(dem, file.path(d, "vat.jpg")))
  expect_equal(c(v$nx, length(v$bands), v$driver), c("1000", "1", "JPEG"))
  st <- rvt_blend(rvt_hillshade(dem), rvt_svf(dem), "multiply", opacity = 0.5)
  s <- img_read(rvt_image(st, file.path(d, "st.webp"), "Viridis"))
  expect_equal(length(s$bands), 4L)
  expect_false(identical(s$bands[[1]], s$bands[[3]]))
})

test_that("relight writes a picture directly", {
  skip_if_offline()
  pt <- c(9.9464, 51.6317)
  crop <- function(p) {
    out <- tempfile(fileext = ".tif")
    gdalraster::translate(p, out, quiet = TRUE, cl_arg = c("-srcwin", "350", "350", "200", "200"))
    out
  }
  rgb <- crop(rvt_data_lgln(pt, "rgb", year = 2025, res = 1))
  dsm <- crop(rvt_data_lgln(pt, "dsm"))
  d <- img_dir()
  on.exit(unlink(c(rgb, dsm, d), recursive = TRUE), add = TRUE)
  for (f in c("r.webp", "r.jpg")) {
    i <- img_read(rvt_relight(rgb, dsm, file.path(d, f)))
    expect_equal(c(i$nx, i$ny), c(200, 200))
  }
  expect_setequal(list.files(d), c("r.webp", "r.jpg"))
})
