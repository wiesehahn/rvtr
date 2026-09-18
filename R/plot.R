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
#' @param path path to a single-band raster, or an `rvt_stack` from
#'   [rvt_blend()], which is rendered to a temporary raster first
#' @param max_dim longest side, in pixels, to render at (default 1000);
#'   downsampling happens on read via GDAL, not by loading the full raster
#'   first
#' @param col a palette name such as `"Grays"` (default), `"Viridis"` or
#'   `"Blue-Red 3"`, or a vector of colours. Names ignore case, spaces and
#'   hyphens, and a unique start is enough. [rvt_palettes()] shows them all.
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
#' p |> rvt_plot(minmax_def = c(-0.05, 0.05))
#'
#' # or cut percentiles, when you don't know the range in advance
#' p |> rvt_plot(minmax_pct_cut = c(2, 98))
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
#' p |> rvt_plot("Blue-Red 3", minmax_def = c(-0.3, 0.3), legend = TRUE)
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
#' * `col` takes a palette name or colours. Use a diverging palette for signed
#'   data, a **qualitative** one for [rvt_geomorphons()] - its values are
#'   class codes, so a ramp implies an order that isn't there:
#'   `rvt_plot(g, "Dark 3")`. [rvt_palettes()] draws every option.
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
#' dem |> rvt_plot()
#' dem |> rvt_plot(legend = TRUE, axes = TRUE, main = basename(dem))
#'
#' # any palette by name; rvt_palettes() shows them all
#' dem |> rvt_svf() |> rvt_plot("Viridis")
#'
#' # signed data: diverging palette, symmetric stretch
#' dem |> rvt_slrm() |> rvt_plot("Blue-Red 3", minmax_def = c(-1, 1), legend = TRUE)
#'
#' # percentile stretch when the range isn't known
#' dem |> rvt_curvature() |> rvt_plot("Blue-Red 3", minmax_pct_cut = c(2, 98))
#'
#' # a blend stack draws directly, with no rvt_render() call
#' dem |> rvt_hillshade() |>
#'   rvt_blend(rvt_svf(dem), "multiply", opacity = 0.25) |>
#'   rvt_plot()
#' @export
rvt_plot <- function(path, col = "Grays", max_dim = 1000, legend = FALSE,
                      axes = FALSE, band = NULL, ...) {
  # a blend stack is a recipe, not a file - render it to a temporary raster so
  # exploring a composite needs no rvt_render() call
  if (inherits(path, "rvt_stack")) path <- rvt_render(path)
  path <- .as_path(path)
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
                           # an RGB composite has no colour map, so a palette
                           # name is only resolved when it will be used
                           col_map_fn = if (length(band) == 1L) .resolve_col(col) else NULL,
                           legend = legend && length(band) == 1L, axes = axes,
                           xlab = "", ylab = "", ...)
  invisible(path)
}
