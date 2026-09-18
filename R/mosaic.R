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
  vrt_path
}
