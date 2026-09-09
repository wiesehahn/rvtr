## Quick-look raster plotting.

#' Quick-look plot of a raster
#'
#' Renders a single-band raster - typically output from [rvt_svf()],
#' [rvt_openness()] or [rvt_vat()] - as an image. A thin wrapper around
#' [gdalraster::plot_raster()], which does the actual work: correct
#' orientation, a colour scale stretched to the data's own range, and NoData
#' handling. This just adds a downsampled read (so a huge raster still
#' renders fast) and a path-in, path-out interface for piping.
#'
#' @param path path to a single-band raster
#' @param max_dim longest side, in pixels, to render at (default 1000);
#'   downsampling happens on read via GDAL, not by loading the full raster
#'   first
#' @param col colour ramp (default: a grey ramp), passed to
#'   [gdalraster::plot_raster()]'s `col_map_fn`
#' @param legend show a colour scale (default `FALSE`, just the image)
#' @param axes show axes with coordinates (default `FALSE`)
#' @param ... passed on to [gdalraster::plot_raster()], e.g. `main = "..."`
#'   for a title
#' @return `path`, invisibly. `path` is the first argument, so this is
#'   pipe-friendly: `path |> rvt_plot()`.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_plot(dem)
#' rvt_plot(dem, legend = TRUE, axes = TRUE, main = basename(dem))
#' @export
rvt_plot <- function(path, max_dim = 1000,
                      col = grDevices::hcl.colors(256, "Grays"),
                      legend = FALSE, axes = FALSE, ...) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  scale <- min(1, max_dim / max(nx, ny))

  gdalraster::plot_raster(ds,
                           xsize = max(1L, round(nx * scale)),
                           ysize = max(1L, round(ny * scale)),
                           col_map_fn = col, legend = legend, axes = axes,
                           xlab = "", ylab = "", ...)
  invisible(path)
}
