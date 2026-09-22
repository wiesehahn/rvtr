# Shared test helpers. testthat sources every helper-*.R before each test file,
# so anything defined here is visible to all of them - and under parallel
# execution, once per subprocess.
#
# Two things live here that are not just convenience:
#
#  - `skip_unless_network()`. The tests that fetch from LGLN and Mapterhorn are
#    about a third of the suite's runtime and were previously gated only by
#    skip_if_offline(), i.e. they always ran. They are opt-in now.
#  - the `fx_*()` fixtures. A hillshade or a sky-view factor of the bundled tile
#    is used as "some finished raster" by a couple of dozen blocks across five
#    files, and each one used to recompute it. These build on first use and
#    cache. **Callers must not unlink() them.**

ref_path <- function(f) system.file("extdata", f, package = "rvtr")
has_ref <- function(f) nzchar(ref_path(f))

dem_path <- function() ref_path("dtm1.tif")
dem <- dem_path()

## The canonical reader: one band as a matrix, NoData as NA. Replaces the
## four near-identical copies that used to sit at the top of test-reference.R,
## test-blend.R and test-image.R.
read_band <- function(p, b = 1L) {
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

# ---------------------------------------------------------------- fixtures --

## Build once, reuse everywhere. Keyed by name; the value is a path that stays
## valid for the rest of the session, so nothing may delete it.
.rv_fx <- new.env(parent = emptyenv())
rv_fx <- function(name, build) {
  if (is.null(.rv_fx[[name]])) .rv_fx[[name]] <- build()
  .rv_fx[[name]]
}

fx_hillshade <- function() rv_fx("hillshade", function() rvt_hillshade(dem_path()))
fx_svf <- function() rv_fx("svf", function() rvt_svf(dem_path(), reach = 10))
fx_slope <- function() rv_fx("slope", function() rvt_slope(dem_path()))
fx_msrm <- function() rv_fx("msrm", function() rvt_msrm(dem_path()))

## Curvature of the bundled tile, keyed by type and radius. Several blocks in
## test-surface.R need the same type twice - "profile" at radius 5 was computed
## three times over - and cross-checking one type against another is the point
## of half of them, so they all come through here.
fx_curvature <- function(type = "profile", radius = NULL) {
  rv_fx(paste0("curv:", type, ":", if (is.null(radius)) "3x3" else radius),
        function() {
          if (is.null(radius)) rvt_curvature(dem_path(), type = type)
          else rvt_curvature(dem_path(), type = type, radius = radius)
        })
}

## The exact full-resolution openness of the 800 px crop at reach 200 - a
## pyramid_px covering the whole reach in one level. It is the reference for
## both "a multi-resolution search matches a full-resolution one" and "a
## coarser pyramid is worse but still close", which used to compute it twice;
## at ~3.3 s it was the most expensive duplicated call in the suite.
fx_opns_exact <- function() {
  rv_fx("opns_exact", function()
    rvt_openness(crop_raster(), reach = 200, pyramid_px = 200))
}

## A cropped copy, for the tests that compare the multi-resolution search
## against an exact full-resolution one. That baseline - pyramid_px equal to
## the whole reach - is what costs, not the pyramid, so the saving comes from
## the raster rather than from weakening the approximation: pyramid_px has to
## stay at its default or the error being measured is a different one.
## The interior margin used with these has to *exceed* the reach, not equal
## it, or the comparison picks up the edge extrapolation the two methods
## deliberately differ on: at 600 px with a 200 px margin the sky-view error
## reads 1.2e-4 against 3.4e-5 well clear of the edge.
crop_raster <- function(src = dem_path(), size = 800L) {
  rv_fx(paste0("crop:", src, ":", size), function() {
    p <- tempfile(fileext = ".tif")
    gdalraster::translate(src, p, quiet = TRUE,
                           cl_arg = c("-srcwin", "0", "0", as.character(size),
                                       as.character(size)))
    p
  })
}

## The bundled tile cut into its four 500x500 quadrants, as four separate
## files. Built once and shared by the mosaic and multi-resolution mosaic
## tests, which used to build their own copies.
quad_rasters <- function() {
  rv_fx("quads", function() {
    vapply(1:4, function(i) {
      p <- tempfile(fileext = ".tif")
      gdalraster::translate(dem_path(), p, quiet = TRUE, cl_arg = c("-srcwin",
        as.character(c(0, 500)[1 + (i - 1) %% 2]),
        as.character(c(0, 500)[1 + (i - 1) %/% 2]), "500", "500"))
      p
    }, "")
  })
}

# ------------------------------------------------------------ alpine tiles --

# Four adjacent 1 km DTM tiles around the Partnachklamm - 494 m of relief,
# against dtm1.tif's 71 m, and a genuine four-file mosaic rather than one
# raster cut into quadrants. Gitignored like the rvt-py references; see
# inst/extdata/README.md for the one-command download.
alp_tiles <- c("dgm1_659_5259.tif", "dgm1_660_5259.tif",
               "dgm1_659_5258.tif", "dgm1_660_5258.tif")
has_alps <- all(vapply(alp_tiles, has_ref, logical(1)))
alps <- function() vapply(alp_tiles, ref_path, "")

# ----------------------------------------------------------------- network --

## Fetching from LGLN and Mapterhorn is roughly a third of the suite's wall
## clock, and far more on a cold cache. Those tests are opt-in: set
## RVTR_TEST_NETWORK=true to run them. Nothing else in the suite needs a
## network connection.
skip_unless_network <- function() {
  if (!isTRUE(as.logical(Sys.getenv("RVTR_TEST_NETWORK", "false"))))
    testthat::skip("network test; set RVTR_TEST_NETWORK=true to run")
  testthat::skip_if_offline()
}
