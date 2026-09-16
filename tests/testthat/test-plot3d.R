# rvt_plot3d(). The picture itself is for looking at, so these check the
# geometry handed to rgl rather than pixels: the grid the DEM is read onto,
# that NoData stays NA, and that a drape is rendered exactly once.

dem_path <- function() system.file("extdata", "dtm1.tif", package = "rvtr")

test_that("arguments are validated before any work", {
  skip_if_not_installed("rgl")
  expect_error(rvt_plot3d(dem_path(), exaggeration = 0), "exaggeration")
  expect_error(rvt_plot3d(dem_path(), exaggeration = c(1, 2)), "exaggeration")
  expect_error(rvt_plot3d(dem_path(), max_dim = 1), "max_dim")
})

test_that("a raster passed as col says so, in words", {
  skip_if_not_installed("rgl")
  # `col` is the second argument in both plot functions, so the slip in this
  # direction is passing the drape positionally. Without the guard it fails
  # inside the palette lookup, which explains nothing.
  hs <- rvt_hillshade(dem_path())
  on.exit(unlink(hs), add = TRUE)
  expect_error(rvt_plot3d(dem_path(), hs), "drape")
  expect_error(rvt_plot3d(dem_path(), rvt_stack(hs)), "drape")
})

test_that("the DEM is read onto a downsampled grid of cell centres", {
  g <- rvtr:::.plot3d_grid(dem_path(), max_dim = 100)
  expect_equal(c(g$nx, g$ny), c(100, 100))
  expect_equal(dim(g$z), c(100, 100))
  expect_length(g$x, 100)
  expect_length(g$y, 100)

  # cell centres sit half a cell inside the bounding box. y ascends, because
  # surface3d() indexes z[i, j] at (x[i], y[j]) and wants y increasing - the
  # raster's north-to-south rows are reversed on read
  cell <- (g$bbox[3] - g$bbox[1]) / 100
  expect_equal(g$x[1], g$bbox[1] + cell / 2)
  expect_equal(g$y[1], g$bbox[2] + cell / 2)
  expect_true(g$y[1] < g$y[100])
})

test_that("z is indexed [x, y] and matches the raster, on a non-square grid", {
  # The invariant that pins the axis order, and it must use a NON-square
  # raster. surface3d() wants z[i, j] at (x[i], y[j]) with y ascending, while
  # read_ds() returns rows north first, west to east - so the naive matrix()
  # fill produces the transpose. On a square grid the transpose has the same
  # dimensions and a test written in the same convention cancels the error on
  # both sides: an earlier version of this test passed while the rendered
  # surface was visibly torn into folds. A rectangle cannot hide it.
  rect <- tempfile(fileext = ".tif")
  gdalraster::translate(dem_path(), rect, quiet = TRUE,
                        cl_arg = c("-srcwin", "0", "0", "800", "400"))
  on.exit(unlink(rect), add = TRUE)

  g <- rvtr:::.plot3d_grid(rect, max_dim = 80)
  expect_equal(c(g$nx, g$ny), c(80, 40))
  expect_equal(dim(g$z), c(80, 40))
  expect_length(g$x, 80)
  expect_length(g$y, 40)
  expect_true(all(diff(g$y) > 0))
  expect_true(all(diff(g$x) > 0))

  # elevations read straight from the raster, not through the grid's own
  # convention, so the comparison is independent
  ds <- methods::new(gdalraster::GDALRaster, rect, read_only = TRUE)
  on.exit(ds$close(), add = TRUE)
  gt <- ds$getGeoTransform()
  at <- function(xi, yj)
    as.numeric(ds$read(1, floor((xi - gt[1]) / gt[2]),
                       floor((yj - gt[4]) / gt[6]), 1, 1, 1, 1))

  probe <- expand.grid(i = c(1L, 2L, 40L, 79L, 80L), j = c(1L, 2L, 20L, 39L, 40L))
  err <- mapply(function(i, j) abs(g$z[i, j] - at(g$x[i], g$y[j])),
                probe$i, probe$j)
  # downsampling averages, so allow a little, but a transposed or flipped grid
  # is out by tens of metres here
  expect_lt(max(err, na.rm = TRUE), 3)
})

test_that("a raster smaller than max_dim is read whole, not upsampled", {
  g <- rvtr:::.plot3d_grid(dem_path(), max_dim = 5000)
  expect_equal(c(g$nx, g$ny), c(1000, 1000))
})

test_that("NoData arrives as NA, so holes become gaps rather than spikes", {
  g <- rvtr:::.plot3d_grid(dem_path(), max_dim = 1000)
  expect_true(anyNA(g$z))
  expect_true(all(g$z > 0, na.rm = TRUE))   # no -9999 or NaN sentinels left
})

test_that("a single-band drape becomes a PNG the size of the read", {
  png <- rvtr:::.plot3d_texture(dem_path(), "Grays", max_dim = 120)
  on.exit(unlink(png))
  ds <- methods::new(gdalraster::GDALRaster, png, read_only = TRUE)
  on.exit(ds$close(), add = TRUE)
  expect_equal(c(ds$getRasterXSize(), ds$getRasterYSize()), c(120, 120))
  # RGB, not a paletted PNG: a colour table makes a poor texture, and
  # grDevices::png() produces one whenever the image has few enough colours
  expect_gte(ds$getRasterCount(), 3L)
  expect_equal(NROW(ds$getColorTable(1)), 0L)
})

test_that("a named palette colours a single-band drape", {
  # this path shipped broken once: .resolve_col() returns the colours
  # themselves, not a ramp function, and calling it as one failed only when a
  # drape was actually coloured - every other test still passed
  for (pal in c("Grays", "Viridis", "Terrain 2")) {
    png <- rvtr:::.plot3d_texture(dem_path(), pal, max_dim = 40)
    expect_true(file.exists(png), label = pal)
    unlink(png)
  }
  expect_error(rvtr:::.plot3d_texture(dem_path(), "nosuchpalette", max_dim = 40),
               "not a palette name")

  # and the same for the no-drape path, which colours the terrain itself
  expect_length(rvtr:::.resolve_col("Grays", 256L), 256L)
})

test_that("the texture keeps the raster's orientation, north at the top", {
  # Checked against the raster itself rather than a rendered frame: OpenGL
  # snapshots proved an unreliable witness here. A hillshade looks plausible
  # whichever way up it is, so the second half uses an asymmetric marker,
  # which no flip can fake.
  n <- 100L
  tex <- rvtr:::.plot3d_texture(dem_path(), "Grays", max_dim = n)
  on.exit(unlink(tex), add = TRUE)

  src <- methods::new(gdalraster::GDALRaster, dem_path(), read_only = TRUE)
  v <- gdalraster::read_ds(src, out_xsize = n, out_ysize = n)
  src$close()
  rng <- range(v, na.rm = TRUE)
  expected <- matrix((v - rng[1]) / diff(rng), n, n, byrow = TRUE)

  d <- methods::new(gdalraster::GDALRaster, tex, read_only = TRUE)
  got <- matrix(d$read(1, 0, 0, n, n, n, n), n, n, byrow = TRUE) / 255
  d$close()

  cmp <- function(a, b) stats::cor(as.numeric(a), as.numeric(b),
                                   use = "complete.obs")
  expect_gt(cmp(got, expected), 0.95)
  expect_lt(cmp(got, expected[n:1, ]), 0.5)      # not flipped north-south
  expect_lt(cmp(got, expected[, n:1]), 0.5)      # not flipped east-west
  expect_lt(cmp(got, t(expected)), 0.5)          # not transposed

  # an unmistakable marker in the north-west of the extent
  mark <- tempfile(fileext = ".tif")
  gdalraster::translate(dem_path(), mark, quiet = TRUE,
                        cl_arg = c("-ot", "Byte", "-scale"))
  ds <- methods::new(gdalraster::GDALRaster, mark, read_only = FALSE)
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  ds$write(1, 0, 0, nx %/% 3, ny %/% 3, rep(255, (nx %/% 3) * (ny %/% 3)))
  ds$close()
  t2 <- rvtr:::.plot3d_texture(mark, "Grays", max_dim = 60)
  on.exit(unlink(c(mark, t2)), add = TRUE)

  d2 <- methods::new(gdalraster::GDALRaster, t2, read_only = TRUE)
  m <- matrix(d2$read(1, 0, 0, 60, 60, 60, 60), 60, 60, byrow = TRUE)
  d2$close()
  quad <- c(NW = mean(m[1:30, 1:30]), NE = mean(m[1:30, 31:60]),
            SW = mean(m[31:60, 1:30]), SE = mean(m[31:60, 31:60]))
  expect_equal(names(which.max(quad)), "NW")
})

test_that("the drape lands the right way up on the surface", {
  skip_if_not_installed("rgl")
  # Nothing here opens an rgl window: the rgl calls are mocked and the
  # arguments rvt_plot3d() hands to surface3d() are inspected directly.
  #
  # What this pins, and why it is the right thing to pin: rgl maps texture
  # t = 1 to the PNG's first row, and .plot3d_texture() writes the north edge
  # as that first row, while the grid's y ascends south to north. So t must
  # ascend with y. That direction was established once by measurement - a
  # north-dark/south-bright drape rendered from a genuine top-down camera gave
  # north 37, south 242 - and CLAUDE.md records it. Flipping t here reverses
  # the drape, which is what this test exists to catch.
  dem <- dem_path()
  hs <- rvt_hillshade(dem)
  on.exit(unlink(hs), add = TRUE)

  args <- NULL
  testthat::with_mocked_bindings(
    rvt_plot3d(dem, drape = hs, max_dim = 60),
    open3d   = function(...) invisible(1L),
    bg3d     = function(...) invisible(NULL),
    aspect3d = function(...) invisible(NULL),
    surface3d = function(...) { args <<- list(...); invisible(1L) },
    .package = "rgl")

  expect_false(is.null(args))
  s <- args$texture_s
  t <- args$texture_t
  expect_equal(dim(s), c(60L, 60L))
  expect_equal(dim(t), c(60L, 60L))

  # s runs west to east across the grid's first index
  expect_equal(s[1, 1], 0)
  expect_equal(s[60, 1], 1)
  expect_true(all(diff(s[, 1]) > 0))

  # t ascends with y: 0 at the south edge, 1 at the north
  expect_equal(t[1, 1], 0)
  expect_equal(t[1, 60], 1)
  expect_true(all(diff(t[1, ]) > 0))

  # and the surface really is the grid, in [x, y] order with y ascending
  g <- rvtr:::.plot3d_grid(dem, max_dim = 60)
  expect_equal(args[[1]], g$x)
  expect_equal(args[[2]], g$y)
  expect_equal(dim(args[[3]]), c(60L, 60L))
})

test_that("a flat drape does not divide by a zero range", {
  flat <- tempfile(fileext = ".tif")
  gdalraster::translate(dem_path(), flat, quiet = TRUE,
                        cl_arg = c("-scale", "0", "1", "5", "5"))
  on.exit(unlink(flat), add = TRUE)
  png <- rvtr:::.plot3d_texture(flat, "Grays", max_dim = 40)
  on.exit(unlink(png), add = TRUE)
  expect_true(file.exists(png))
})

test_that("a blend stack drape is rendered once", {
  dem <- dem_path()
  hs <- rvt_hillshade(dem)
  on.exit(unlink(hs), add = TRUE)
  stack <- rvt_stack(hs)

  rendered <- 0L
  png <- testthat::with_mocked_bindings(
    rvtr:::.plot3d_texture(stack, "Grays", max_dim = 80),
    rvt_render = function(...) { rendered <<- rendered + 1L; hs },
    .package = "rvtr")
  on.exit(unlink(png), add = TRUE)
  expect_equal(rendered, 1L)
})

test_that("a three-band drape keeps its own colours", {
  dem <- dem_path()
  hs <- rvt_hillshade(dem)
  rgb3 <- tempfile(fileext = ".tif")
  gdalraster::translate(hs, rgb3, quiet = TRUE,
                        cl_arg = c("-b", "1", "-b", "1", "-b", "1",
                                   "-ot", "Byte", "-scale"))
  on.exit(unlink(c(hs, rgb3)), add = TRUE)

  png <- rvtr:::.plot3d_texture(rgb3, "Grays", max_dim = 60)
  on.exit(unlink(png), add = TRUE)
  ds <- methods::new(gdalraster::GDALRaster, png, read_only = TRUE)
  on.exit(ds$close(), add = TRUE)
  expect_equal(c(ds$getRasterXSize(), ds$getRasterYSize()), c(60, 60))
})
