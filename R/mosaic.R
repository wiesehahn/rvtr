## Mosaic handling.
##
## A horizon search near a tile edge needs terrain from the neighbouring tile,
## so processing a mosaic file-by-file would leave a seam at every boundary.
## Wrapping the whole collection in a GDAL virtual raster removes the problem
## entirely: the VRT is a single addressable dataset, windowed reads span the
## source files transparently, and the tiling engine is then free to cut tiles
## wherever it likes without regard to the original file layout.

#' Present one or many rasters as a single dataset
#'
#' Passes a single path straight through. For several paths, builds a GDAL
#' virtual raster (VRT) over them and returns its path - a small XML file, no
#' pixel data is copied. All functions in this package accept either form, so
#' calling this directly is only necessary if you want to keep and reuse the
#' VRT.
#'
#' @param dem a path, or a character vector of paths to mosaic. A directory,
#'   or a vector of length 1 ending in `.vrt`, is returned unchanged.
#' @param vrt_path where to write the VRT (default: a temp file)
#' @param cl_arg extra arguments passed to `gdalbuildvrt`, e.g.
#'   `c("-resolution", "highest")`
#' @return path to a single dataset usable by [rvt_svf()], [rvt_openness()],
#'   [rvt_vat()]. Pipe-friendly: `tiles |> rvt_mosaic() |> rvt_vat()`.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' dem |> rvt_mosaic()  # single path passes through unchanged
#' @export
rvt_mosaic <- function(dem, vrt_path = fs::file_temp(ext = "vrt"), cl_arg = NULL) {
  dem <- .as_path(dem)
  vrt_path <- .out_path(vrt_path)
  if (length(dem) == 1L) return(dem)
  missing <- dem[!fs::file_exists(dem)]
  if (length(missing))
    stop("input raster(s) not found: ", paste(utils::head(missing, 3), collapse = ", "),
         call. = FALSE)
  gdalraster::buildVRT(vrt_filename = vrt_path, input_rasters = dem,
                        cl_arg = cl_arg, quiet = TRUE)
  # gdalbuildvrt silently leaves out any file whose coordinate system differs
  # from the first one's, so a mosaic could come back with holes and no
  # warning. Check, and only then look closer.
  if (.vrt_n_sources(vrt_path) < length(unique(dem)))
    .mosaic_mixed_crs(dem, vrt_path, cl_arg)
  vrt_path
}

.vrt_n_sources <- function(vrt) {
  x <- paste(readLines(vrt, warn = FALSE), collapse = "\n")
  length(unique(regmatches(x, gregexpr("<SourceFilename[^>]*>[^<]*</SourceFilename>",
                                       x))[[1]]))
}

## Files that differ only in a vertical CRS - LGLN's autumn-2025 bDOM carries
## ETRS89 / UTM 32N + DHHN2016 beside plain EPSG:25832 tiles - are the same
## grid and are mosaicked under the common horizontal CRS. Anything else is a
## genuinely different projection, and an error rather than a silent hole.
.mosaic_mixed_crs <- function(dem, vrt_path, cl_arg) {
  srs <- vapply(dem, function(p) {
    ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
    on.exit(ds$close())
    ds$getProjectionRef()
  }, character(1))
  h <- vapply(srs, .srs_horizontal, character(1))
  same <- vapply(h, function(x) {
    if (!nzchar(x) || !nzchar(h[1])) return(nzchar(x) == nzchar(h[1]))
    gdalraster::srs_is_same(x, h[1])
  }, logical(1))
  if (all(same)) {
    gdalraster::buildVRT(vrt_filename = vrt_path, input_rasters = dem, quiet = TRUE,
      cl_arg = c("-allow_projection_difference",
                 if (nzchar(h[1])) c("-a_srs", h[1]), cl_arg))
    if (.vrt_n_sources(vrt_path) == length(unique(dem))) return(invisible(vrt_path))
  }
  odd <- if (all(same)) dem else dem[!same]
  stop(sprintf(paste("Can't mosaic these rasters: %s %s a different %s from the",
                     "first one (%s). Bring them onto one grid first, e.g. with",
                     "gdalwarp."),
               paste(utils::head(basename(odd), 3), collapse = ", "),
               if (length(odd) == 1L) "has" else "have",
               if (all(same)) "band layout or data type" else "coordinate system",
               basename(dem[1])), call. = FALSE)
}
