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
#' @param band which band to draw, for multi-band results such as
#'   [rvt_multi_hillshade()]; or three band numbers to draw as red, green and
#'   blue. Defaults to band 1, except for a 3-band raster, which is drawn as
#'   RGB (which is what [rvt_mstp()] wants).
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
                      legend = FALSE, axes = FALSE, band = NULL, ...) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  scale <- min(1, max_dim / max(nx, ny))
  out_nx <- max(1L, round(nx * scale))
  out_ny <- max(1L, round(ny * scale))

  if (is.null(band)) band <- if (ds$getRasterCount() == 3L) 1:3 else 1L

  # read_ds() tags the values with the georeferencing plot_raster() needs, so
  # picking a band here still gives the same picture as handing over the whole
  # dataset would.
  v <- gdalraster::read_ds(ds, bands = band, out_xsize = out_nx,
                            out_ysize = out_ny)

  gdalraster::plot_raster(v, xsize = out_nx, ysize = out_ny,
                           nbands = length(band),
                           col_map_fn = if (length(band) == 1L) col else NULL,
                           legend = legend && length(band) == 1L, axes = axes,
                           xlab = "", ylab = "", ...)
  invisible(path)
}
