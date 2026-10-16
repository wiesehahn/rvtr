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
#' @param range `c(lo, hi)` to stretch the colours over, or `NULL` (default)
#'   for the data's own minimum and maximum. Three bands share one range; see
#'   **Controlling the stretch**.
#' @param pct `c(low, high)` percentages to cut from each tail instead, e.g.
#'   `c(2, 98)`, measured with [rvt_range()]. Not with `range`.
#' @param stretch `"auto"` (default) uses `"tanh"` for signed data and
#'   `"linear"` for everything else; `"linear"` clips everything outside the
#'   range to the end colours; `"tanh"` compresses it softly instead - see
#'   **Soft stretch**
#' @param contrast with `stretch = "tanh"`, how steep the curve is (default
#'   1.5): higher brings out more mid-range texture and pushes the tails
#'   further towards the end colours
#' @param ... passed on to [gdalraster::plot_raster()] - see **Controlling the
#'   stretch** and **Other display options** below
#' @return `path`, invisibly. `path` is the first argument, so this is
#'   pipe-friendly: `path |> rvt_plot()`.
#'
#' @section Controlling the stretch:
#' Signed data gets a soft stretch by default (see **Soft stretch**).
#' Everything else is stretched linearly across the data's full range, which is
#' often the wrong choice: a handful of extreme cells then compress everything
#' interesting into the middle of the palette. Two arguments override it, and
#' they are spelled as in [rvt_image()], so a stretch chosen on screen carries
#' over to the written picture unchanged:
#'
#' ```r
#' # fixed values - anything outside is clipped to the end colours
#' p |> rvt_plot(range = c(-0.05, 0.05))
#'
#' # or cut percentiles, when you don't know the range in advance
#' p |> rvt_plot(pct = c(2, 98))
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
#' p |> rvt_plot("Blue-Red 3", range = c(-0.3, 0.3), legend = TRUE)
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
#' Three bands drawn as RGB share one `range`, as in [rvt_image()], so the
#' colour balance is not shifted band by band: `range = c(0, 255)` for an
#' orthophoto. [gdalraster::plot_raster()]'s own `minmax_def` and
#' `minmax_pct_cut` still work through `...` when you do want a separate
#' stretch per band, but they cannot be combined with `range` or `pct`.
#'
#' @section Soft stretch:
#' A linear stretch clips: every value beyond `range` gets the end colour, so
#' the inside of the deepest gully or the top of the highest bank reads as one
#' flat black or white. `stretch = "tanh"` bends the colour scale into an S
#' instead. `range` (or `pct`) still sets where the contrast is concentrated,
#' but its ends land at 5% and 95% of the palette rather than 0% and 100%, and
#' values beyond fade gradually into the end colours. Detail survives in the
#' tails, and the middle gains a little contrast.
#'
#' ```r
#' p |> rvt_plot(pct = c(2, 98), stretch = "tanh")
#' p |> rvt_plot(pct = c(2, 98), stretch = "tanh", contrast = 2.5)  # harder
#' ```
#'
#' It suits signed results - [rvt_msrm()], [rvt_slrm()], [rvt_tpi()],
#' [rvt_dev()], [rvt_log()], curvature. On a lopsided one such as [rvt_svf()],
#' where most cells sit close to 1, the curve is centred in the wrong place and
#' mainly washes out the open ground. So **the default, `stretch = "auto"`,
#' decides from the data**: a single band whose 2-98% range crosses zero, with
#' neither tail more than four times the other, is drawn with `"tanh"`, and
#' everything else - positive data, class codes, RGB - linearly. Name the
#' stretch to override it.
#'
#' A tanh range that crosses zero is **centred on zero**, spreading as far as
#' the larger tail, so a diverging palette's neutral colour sits on zero and
#' convex mirrors concave even when the two tails differ. Without `range` or
#' `pct`, tanh takes the 2-98% range rather than the full minimum and maximum,
#' which would centre the curve on a few outliers.
#'
#' With a legend, the scale stays in data units and shows the curve.
#' [rvt_image()] takes the same arguments and decides the same way, so a
#' stretch settled on screen is the one you write.
#'
#' The idea is the tone map of Leland Brown's texture shading, as implemented
#' in `texture-shade`.
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
#' dem |> rvt_slrm() |> rvt_plot("Blue-Red 3", range = c(-1, 1), legend = TRUE)
#'
#' # percentile stretch when the range isn't known
#' dem |> rvt_curvature() |> rvt_plot("Blue-Red 3", pct = c(2, 98))
#'
#' # signed data gets soft tails by default; ask for linear to clip instead
#' dem |> rvt_msrm() |> rvt_plot("Blue-Red 3")
#' dem |> rvt_msrm() |> rvt_plot("Blue-Red 3", pct = c(2, 98), stretch = "linear")
#'
#' # a blend stack draws directly, with no rvt_render() call
#' dem |> rvt_hillshade() |>
#'   rvt_blend(rvt_svf(dem), "multiply", opacity = 0.25) |>
#'   rvt_plot()
#' @export
rvt_plot <- function(path, col = "Grays", max_dim = 1000, legend = FALSE,
                      axes = FALSE, band = NULL, range = NULL, pct = NULL,
                      stretch = "auto", contrast = 1.5, ...) {
  stretch <- .check_stretch(stretch, contrast)
  # `range`/`pct` are rvt_image()'s spelling of plot_raster()'s minmax_def and
  # minmax_pct_cut, resolved here to one range shared by every band drawn -
  # per-band stretching would shift an RGB image's colour balance. Passing
  # gdalraster's own arguments through ... still gives per-band control.
  dots <- list(...)
  gd_stretch <- any(c("minmax_def", "minmax_pct_cut") %in% names(dots))
  if (stretch == "auto" && gd_stretch) stretch <- "linear"
  if ((!is.null(range) || !is.null(pct) || stretch == "tanh") && gd_stretch)
    stop(paste("Give either `range`/`pct`/`stretch` or `minmax_def`/`minmax_pct_cut`,",
               "not both."), call. = FALSE)
  if (!is.null(range) && !is.null(pct))
    stop("Give either `range` or `pct`, not both.", call. = FALSE)
  if (!is.null(range) && (length(range) != 2L || !all(is.finite(range)) ||
                           range[1] >= range[2]))
    stop("`range` must be `c(lo, hi)` with lo < hi, or NULL", call. = FALSE)

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

  # one range across all bands drawn, measured the same way rvt_image() does
  if (!is.null(pct)) range <- .pct_range(path, band, pct)
  # an RGB composite has no colour map, so a palette name is only resolved
  # when it will be used
  col_map <- if (length(band) == 1L) .resolve_col(col) else NULL

  q98 <- NULL
  if (stretch == "auto") {
    a <- .auto_stretch(path, band)
    stretch <- a$stretch; q98 <- a$q
  }

  if (stretch == "tanh") {
    # plot_raster() hands col_map_fn values already stretched to 0-1 over
    # minmax_def, and draws the legend through the same function. So the
    # stretch spans the data's full extent, and the curve goes into the colour
    # function: the soft tails stay visible and the legend stays in data units.
    full <- vapply(band, function(b) rvt_range(path, band = b), numeric(2))
    full <- c(min(full[1, ]), max(full[2, ]))
    # without a range, 2-98%: centred on the full extent the curve would sit
    # on the outliers and leave the bulk of the data a flat mid-tone
    if (is.null(range)) range <- if (!is.null(q98)) q98 else .pct_range(path, band)
    r <- .tanh_range(range)
    dots$minmax_def <- c(rep(full[1], length(band)), rep(full[2], length(band)))
    soft <- function(a) .norm_tanh(full[1] + a * (full[2] - full[1]),
                                   r[1], r[2], contrast)
    col_map <- .tanh_col_map(soft, col_map, length(band))
  } else if (!is.null(range)) {
    dots$minmax_def <- c(rep(range[1], length(band)), rep(range[2], length(band)))
  }

  # read_ds() tags the values with the georeferencing plot_raster() needs, so
  # picking a band here still gives the same picture as handing over the whole
  # dataset would.
  v <- gdalraster::read_ds(ds, bands = band, out_xsize = out_nx,
                            out_ysize = out_ny)

  do.call(gdalraster::plot_raster,
          c(list(v, xsize = out_nx, ysize = out_ny,
                 nbands = length(band),
                 col_map_fn = col_map,
                 legend = legend && length(band) == 1L, axes = axes,
                 xlab = "", ylab = ""),
            dots))
  invisible(path)
}

## A col_map_fn applying the tanh curve before the colours: `soft` takes
## plot_raster()'s 0-1 values and returns the 0-1 position on the palette.
## NoData stays NA, which plot_raster() then draws in its own na_col.
.tanh_col_map <- function(soft, cols, nbands) {
  if (nbands == 1L) {
    if (is.function(cols)) return(function(a) cols(soft(a)))
    ramp <- grDevices::colorRamp(cols)
    return(function(a) {
      out <- rep(NA_character_, length(a))
      ok <- !is.na(a)
      out[ok] <- grDevices::rgb(ramp(soft(a[ok])), maxColorValue = 255)
      out
    })
  }
  function(r, g, b) {
    out <- rep(NA_character_, length(r))
    ok <- !(is.na(r) | is.na(g) | is.na(b))
    out[ok] <- grDevices::rgb(soft(r[ok]), soft(g[ok]), soft(b[ok]))
    out
  }
}
