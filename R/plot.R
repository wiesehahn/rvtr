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
#'   [rvt_mstp()]; or three band numbers to draw as red, green and blue.
#'   Defaults to band 1, except for a 3-band raster, which is drawn as RGB
#'   (which is what [rvt_mstp()] wants).
#' @param ... passed on to [gdalraster::plot_raster()] - see **Controlling the
#'   stretch** and **Other display options** below
#' @return `path`, invisibly. `path` is the first argument, so this is
#'   pipe-friendly: `path |> rvt_plot()`.
#'
#' @section Controlling the stretch:
#' By default the colours are stretched across the data's full range, which is
#' often the wrong choice: a handful of extreme cells then compress everything
#' interesting into the middle of the palette. Two arguments override it,
#' both passed through to [gdalraster::plot_raster()]:
#'
#' ```r
#' # fixed values - anything outside is clipped to the end colours
#' rvt_plot(p, minmax_def = c(-0.05, 0.05))
#'
#' # or cut percentiles, when you don't know the range in advance
#' rvt_plot(p, minmax_pct_cut = c(2, 98))
#' ```
#'
#' The difference is easy to underestimate. Minimal curvature on the bundled
#' sample tile spans -1.41 to 0.87, but its 2nd-98th percentile range is only
#' -0.20 to 0.12 - so the default stretch spends most of the palette on
#' outliers.
#'
#' For any **signed** result - curvature, [rvt_log()], [rvt_slrm()],
#' [rvt_msrm()], [rvt_dev()] - keep the range **symmetric about zero** and use
#' a diverging palette, or zero will not land on the palette's neutral centre
#' and convex will stop reading as the opposite of concave:
#'
#' ```r
#' div <- grDevices::hcl.colors(256, "Blue-Red 3")
#' p |> rvt_plot(col = div, minmax_def = c(-0.3, 0.3), legend = TRUE)
#' ```
#'
#' Turning the `legend` on while choosing a stretch is worth the space - it is
#' the only way to see what the range you picked actually was. Note the legend
#' is drawn for single-band images only; an RGB composite has no one scale to
#' show.
#'
#' Which numbers to use is metric-dependent, and each metric's own help page
#' carries the stretches recommended in the literature - see the
#' *Recommended settings* section of [rvt_svf()], [rvt_openness()] and the
#' rest.
#'
#' @section Other display options:
#' Anything [gdalraster::plot_raster()] accepts works here. The ones worth
#' knowing:
#'
#' * `main = "..."` for a title; `axes = TRUE` adds coordinate axes.
#' * `col` takes any colour ramp. Use a diverging one for signed data, a
#'   **discrete** one for [rvt_geomorphons()] - its values are class codes, so
#'   a continuous ramp implies an order that isn't there:
#'   `rvt_plot(g, col = grDevices::hcl.colors(10, "Spectral"))`.
#' * `band = i` picks one band out of a multi-band result; a 3-band one such
#'   as [rvt_mstp()] is drawn as RGB automatically.
#' * `max_dim` caps the rendered size. Downsampling happens during the read,
#'   so previewing a huge raster is cheap - and cheaper still on the
#'   Cloud-Optimized GeoTIFFs this package writes, which carry overviews.
#'
#' One caveat when arranging several plots side by side: `legend = TRUE`
#' splits the device with [graphics::layout()], which overrides any
#' `par(mfrow = ...)` you have set, so every panel after the first is lost.
#' Leave the legend off for multi-panel figures, or draw each panel to its own
#' device.
#'
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_plot(dem)
#' rvt_plot(dem, legend = TRUE, axes = TRUE, main = basename(dem))
#'
#' # signed data: diverging palette, symmetric stretch
#' div <- grDevices::hcl.colors(256, "Blue-Red 3")
#' dem |> rvt_slrm() |> rvt_plot(col = div, minmax_def = c(-1, 1), legend = TRUE)
#'
#' # percentile stretch when the range isn't known
#' dem |> rvt_curvature() |> rvt_plot(col = div, minmax_pct_cut = c(2, 98))
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
