## Layer blending: compose finished rasters into a display image.
##
## Deliberately operates on *rasters that already exist*, not fused into the
## metric pass the way rvt_vat() is. Fusing is faster in one shot, but makes
## every opacity tweak recompute the horizon searches - fatal for the
## iterate-until-it-looks-right workflow this exists for. Metrics are computed
## once to COGs; a blend then costs seconds.
##
## rvt_vat() is *not* built on this. It is validated against rvt-py down to a
## deliberate bug flag, and that is worth more than sharing code. The test
## suite asserts this engine reproduces .vat_stack() instead, which pins the
## engine to an already-validated result without touching it.

## RVT's normalize_lin(): linear cut-off then stretch to 0-1. Shared with
## rvt_vat() (R/vat.R), which is where these three used to live.
.norm_lin <- function(x, lo, hi) {
  x[x < lo] <- lo
  x[x > hi] <- hi
  if (hi == lo) return(x * 0)
  (x - lo) / (hi - lo)
}

.blend_overlay <- function(active, background) {
  out <- background
  hi <- !is.na(background) & background > 0.5
  lo <- !is.na(background) & background <= 0.5
  out[hi] <- 1 - (1 - 2 * (background[hi] - 0.5)) * (1 - active[hi])
  out[lo] <- 2 * background[lo] * active[lo]
  out
}

.apply_opacity <- function(active, background, opacity) {
  active * opacity + background * (1 - opacity)
}

## The thirteen QGIS layer blending modes, by their QGIS names so the
## vocabulary transfers. `a` is the layer being laid on top, `b` what is
## already beneath it. Both are 0-1; results are clamped where a mode can
## leave the range.
##
## Asymmetric modes (overlay, dodge, burn, hard_light, subtract, difference is
## not) are why argument order is part of the public contract: the piped
## raster is `b`, the argument is `a`.
##
## `luminosity`, which RVT's VAT names for its slope layer, reduces to
## `normal` on single-band data - there is nothing to implement.
.clamp01 <- function(x) {
  x[!is.na(x) & x < 0] <- 0
  x[!is.na(x) & x > 1] <- 1
  x
}

.blend_modes <- list(
  normal     = function(a, b) a,
  lighten    = function(a, b) pmax(a, b),
  darken     = function(a, b) pmin(a, b),
  multiply   = function(a, b) a * b,
  screen     = function(a, b) 1 - (1 - a) * (1 - b),
  addition   = function(a, b) .clamp01(a + b),
  subtract   = function(a, b) .clamp01(b - a),
  difference = function(a, b) abs(a - b),
  overlay    = function(a, b) .blend_overlay(a, b),
  # hard light is overlay with the two swapped
  hard_light = function(a, b) .blend_overlay(b, a),
  dodge      = function(a, b) {
    out <- ifelse(a >= 1, 1, b / (1 - a))
    .clamp01(out)
  },
  burn       = function(a, b) {
    out <- ifelse(a <= 0, 0, 1 - (1 - b) / a)
    .clamp01(out)
  },
  # W3C / Photoshop soft light
  soft_light = function(a, b) {
    d <- ifelse(b <= 0.25, ((16 * b - 12) * b + 4) * b, sqrt(b))
    .clamp01(ifelse(a <= 0.5,
                    b - (1 - 2 * a) * b * (1 - b),
                    b + (2 * a - 1) * (d - b)))
  }
)

#' Blend modes
#'
#' The thirteen layer blending modes [rvt_blend()] understands, named as they
#' are in QGIS so the vocabulary carries over.
#'
#' @format Character vector of length 13.
#' @seealso [rvt_blend()]
#' @examples
#' rvt_blend_modes
#' @export
rvt_blend_modes <- c("normal", "lighten", "darken", "multiply", "screen",
                      "addition", "subtract", "difference", "overlay",
                      "hard_light", "dodge", "burn", "soft_light")

## A raster's grid: size, cell size, extent and CRS.
.grid <- function(path) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  res <- ds$res()
  list(nx = ds$getRasterXSize(), ny = ds$getRasterYSize(),
       xres = res[1], yres = res[2], bbox = as.numeric(ds$bbox()),
       srs = ds$getProjection())
}

## Two grids can be combined when they cover the same extent in the same CRS
## and one is an exact whole-factor subdivision of the other - 1 m over 0.2 m,
## not 1 m over 0.3 m. Returns the finer grid, which is the one layers are
## brought onto. Checked up front rather than discovered as garbage output;
## comparing extents matters as much as sizes, since two equally sized rasters
## shifted against each other would otherwise blend misaligned.
.grid_join <- function(a, b, what, base = "the stack") {
  fmt <- function(g) sprintf("%d x %d at %g m", g$nx, g$ny, g$xres)
  if (nzchar(a$srs) && nzchar(b$srs) && !gdalraster::srs_is_same(a$srs, b$srs))
    stop(sprintf("`%s` is in a different coordinate system from %s.", what, base),
         call. = FALSE)
  whole <- function(r) r >= 1 && abs(r - round(r)) < 1e-6
  rx <- a$xres / b$xres; ry <- a$yres / b$yres
  finer <- if (whole(rx) && whole(ry)) b
           else if (whole(1 / rx) && whole(1 / ry)) a
           else NULL
  if (is.null(finer))
    stop(sprintf(paste0(
      "`%s` is %s, but %s is %s.\nResolutions can differ, but only by a whole ",
      "factor (e.g. 1 m and 0.2 m) - see `rvt_resample()`."),
      what, fmt(b), base, fmt(a)), call. = FALSE)
  tol <- min(a$xres, b$xres) * 1e-3
  if (any(abs(a$bbox - b$bbox) > tol))
    stop(sprintf(paste0(
      "`%s` covers %s, but %s covers %s.\nLayers must cover the same extent; ",
      "fetch or crop them for the same area."),
      what, paste(format(b$bbox, nsmall = 1), collapse = ", "), base,
      paste(format(a$bbox, nsmall = 1), collapse = ", ")), call. = FALSE)
  finer
}

## Bring a coarser raster onto grid `g` (same extent, whole-factor finer).
## Written to a real file, never left as a lazy VRT, so no read downstream can
## depend on tile_size. Bilinear: smooth, and unlike cubic it cannot overshoot,
## which on a 0/1 shadow layer would ring around every shadow edge. Nearest
## (cell replication) was tried and judged worse on screen, despite keeping
## edges exact. -ovr NONE keeps a source COG's overviews out of it.
.upsample <- function(path, g) {
  out <- fs::file_temp(ext = "tif")
  gdalraster::translate(path, out, quiet = TRUE,
    cl_arg = c("-outsize", g$nx, g$ny, "-r", "bilinear", "-ovr", "NONE",
               "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE",
               "-co", "BIGTIFF=IF_SAFER"))
  out
}

.on_grid <- function(path, g) {
  h <- .grid(path)
  h$nx == g$nx && h$ny == g$ny
}

.nbands <- function(path) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  ds$getRasterCount()
}

#' Stretch range for a layer
#'
#' The `c(lo, hi)` a layer is stretched over before blending, read from the
#' raster itself. Evaluate it when you build the stack, so the numbers are
#' fixed before any tile is touched.
#'
#' That is the whole point: a range derived per tile would give the same ground
#' different values depending on which tile it landed in, which is the one
#' thing this package will not do. Passing a frozen pair sidesteps the question
#' entirely.
#'
#' @param path raster to measure
#' @param pct optional `c(low, high)` percentages to cut from each tail, e.g.
#'   `c(2, 98)`. `NULL` (default) uses the full minimum and maximum.
#' @param band band to measure (default 1)
#' @param max_dim when `pct` is given, the longest side of the decimated read
#'   the percentiles are taken from (default 1000)
#' @return `c(lo, hi)`
#' @seealso [rvt_blend()]
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' svf <- rvt_svf(dem)
#' rvt_range(svf)
#' rvt_range(svf, pct = c(2, 98))
#' @export
rvt_range <- function(path, pct = NULL, band = 1L, max_dim = 1000L) {
  path <- .as_path(path)
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  if (is.null(pct)) {
    s <- ds$getStatistics(band = band, approx_ok = TRUE, force = TRUE)
    return(c(s[1], s[2]))
  }
  if (length(pct) != 2L || any(pct < 0) || any(pct > 100) || pct[1] >= pct[2])
    stop("`pct` must be two increasing percentages in 0-100", call. = FALSE)
  # a decimated read, the same trick rvt_plot() uses: plenty for a percentile
  # and deterministic, so the range does not drift between runs
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  f <- max(1, ceiling(max(nx, ny) / max_dim))
  ox <- max(1L, as.integer(nx %/% f)); oy <- max(1L, as.integer(ny %/% f))
  v <- as.numeric(ds$read(band, 0, 0, nx, ny, ox, oy))
  nd <- ds$getNoDataValue(band)
  if (!is.na(nd)) v[v == nd] <- NA
  unname(stats::quantile(v, pct / 100, na.rm = TRUE))
}

.layer <- function(path, mode, opacity, range, invert) {
  path <- .as_path(path)                     # also accepts a terra SpatRaster
  mode <- match.arg(mode, rvt_blend_modes)
  if (!is.numeric(opacity) || length(opacity) != 1L ||
      is.na(opacity) || opacity < 0 || opacity > 1)
    stop("`opacity` must be a single number in 0-1", call. = FALSE)
  if (!is.null(range) && (length(range) != 2L || !all(is.finite(range)) ||
                           range[1] >= range[2]))
    stop("`range` must be `c(lo, hi)` with lo < hi, or NULL", call. = FALSE)
  list(path = path, mode = mode, opacity = opacity,
       range = range, invert = isTRUE(invert), nbands = .nbands(path))
}

#' Start a blend stack
#'
#' The background layer a stack is built on. Only needed when the background
#' itself needs a `range` or `invert`; otherwise pipe the raster straight into
#' [rvt_blend()], which starts a stack for you.
#'
#' @param base path to the background raster
#' @param range `c(lo, hi)` to stretch over, or `NULL` (default) for the
#'   raster's own minimum and maximum. See [rvt_range()].
#' @param invert flip the layer after stretching, so high reads dark
#' @return an `rvt_stack`
#' @seealso [rvt_blend()] to add layers, [rvt_render()] to write the result
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_stack(rvt_hillshade(dem), range = c(0, 1))
#' @export
rvt_stack <- function(base, range = NULL, invert = FALSE) {
  if (inherits(base, "rvt_stack")) return(base)
  l <- .layer(base, "normal", 1, range, invert)
  s <- list(layers = list(l), grid = .grid(l$path))
  class(s) <- "rvt_stack"
  s
}

#' Blend a layer over a raster
#'
#' Lays one raster over another with a blend mode, building up a display image
#' the way layers stack in QGIS. Pipe rasters through it and render once at the
#' end.
#'
#' ```r
#' hillshade |>
#'   rvt_blend(openness, "overlay",  opacity = 0.5) |>
#'   rvt_blend(svf,      "multiply", opacity = 0.25) |>
#'   rvt_render("relief.tif")
#' ```
#'
#' **Order is bottom to top.** The raster piped in is the *background*; the one
#' passed as an argument goes *over* it. That matters, because `overlay`,
#' `dodge`, `burn`, `hard_light` and `subtract` all treat their two inputs
#' differently. It reads the way a QGIS layer panel does, upwards.
#'
#' **Nothing is computed until you render.** `rvt_blend()` returns a recipe,
#' not a raster, so a chain of layers costs one pass over the data instead of
#' one per layer. [rvt_render()] writes the file; [rvt_plot()] will draw a
#' stack directly, so exploring needs no render call at all.
#'
#' @section Choosing a mode:
#' Three of the thirteen modes cover most composites. `multiply` only ever
#' darkens, so it is what you want for deepening enclosed ground with a
#' sky-view factor. `screen` only ever lightens. `overlay` does both at once -
#' darkening what is already dark and lightening what is already light - which
#' adds contrast without washing the image out, and is why [rvt_vat()] uses it
#' for its openness layer.
#'
#' The rest are situational: `darken`/`lighten` take the extreme of the two,
#' `difference` shows where two visualizations disagree, `soft_light` is a
#' gentler `overlay`, and `dodge`/`burn` are aggressive enough to clip. All
#' thirteen are in [rvt_blend_modes].
#'
#' @section Ranges:
#' Every layer is stretched to 0-1 before blending. Leave `range` alone and the
#' raster's own minimum and maximum are used, measured once before any tile is
#' read. Set it to concentrate contrast where it matters - `c(0.7, 1)` on a
#' sky-view factor spreads out the band where ordinary terrain actually lives.
#' [rvt_range()] reads a percentile stretch off the raster.
#'
#' @section Different resolutions:
#' Layers may differ in resolution by a whole factor - a 0.2 m orthophoto under
#' a 1 m hillshade, say - as long as they cover the same extent in the same
#' coordinate system. The result is rendered at the finest resolution in the
#' stack, and coarser layers are upsampled onto it with bilinear interpolation
#' when rendering. Rasters fetched for one place with [rvt_data_lgln()] or
#' [rvt_data_mapterhorn()] line up this way. Anything else stops with an error
#' naming the mismatch, rather than blending misaligned.
#'
#' A three-band layer is stretched over a single range spanning all its bands,
#' never one range per band: stretching each channel to its own extremes would
#' shift the colour balance. Stretch the channels yourself first if that is
#' what you want.
#'
#' @param x a raster path, or an `rvt_stack` to add to
#' @param layer path to the raster to lay over `x`
#' @param mode one of [rvt_blend_modes]; `"normal"` (default) just replaces,
#'   which with `opacity` is a plain fade between the two
#' @param opacity 0-1, applied after the mode (default 1)
#' @param range `c(lo, hi)` to stretch `layer` over, or `NULL` (default) for
#'   its own minimum and maximum
#' @param invert flip the layer after stretching, so high reads dark. Slope
#'   and negative openness are usually shown this way.
#' @return an `rvt_stack`
#' @seealso [rvt_render()], [rvt_range()], [rvt_blend_modes]; [rvt_vat()] for a
#'   ready-made archaeological composite.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' hs <- rvt_hillshade(dem)
#' svf <- rvt_svf(dem)
#' rvt_blend(hs, svf, "multiply", opacity = 0.25)
#' @export
rvt_blend <- function(x, layer, mode = "normal", opacity = 1,
                       range = NULL, invert = FALSE) {
  s <- rvt_stack(x)
  l <- .layer(layer, mode, opacity, range, invert)
  s$grid <- .grid_join(s$grid, .grid(l$path), deparse(substitute(layer)))
  s$layers <- c(s$layers, list(l))
  s
}

#' @export
print.rvt_stack <- function(x, ...) {
  cat(sprintf("<rvt_stack> %d layer%s, bottom to top\n",
              length(x$layers), if (length(x$layers) == 1L) "" else "s"))
  for (i in seq_along(x$layers)) {
    l <- x$layers[[i]]
    cat(sprintf("  %d. %-28s %-11s opacity %-5.2f %s%s\n", i,
                basename(l$path), if (i == 1L) "(background)" else l$mode,
                l$opacity,
                if (is.null(l$range)) "range auto"
                else sprintf("range %g-%g", l$range[1], l$range[2]),
                if (l$invert) ", inverted" else ""))
  }
  g <- x$grid
  cat(sprintf("Renders at %d x %d, %g m. Not a raster yet - call rvt_render() to write it.\n",
              g$nx, g$ny, g$xres))
  invisible(x)
}

## One tile: normalise each layer, blend it over what is beneath, apply
## opacity. `aux` entries arrive in layer order, one per band per layer.
.blend_tile <- function(tile, layers, nb, ctx) {
  # the background, broadcast to the output band count
  base <- layers[[1]]
  bg <- vector("list", nb)
  for (k in seq_len(nb)) {
    m <- if (base$nbands == 1L && k > 1L) ctx$aux[[base$aux[1]]]$m
         else ctx$aux[[base$aux[min(k, base$nbands)]]]$m
    v <- .norm_lin(m, base$range[1], base$range[2])
    bg[[k]] <- if (base$invert) 1 - v else v
  }

  for (i in seq_along(layers)[-1]) {
    l <- layers[[i]]
    f <- .blend_modes[[l$mode]]
    for (k in seq_len(nb)) {
      m <- ctx$aux[[l$aux[min(k, l$nbands)]]]$m
      a <- .norm_lin(m, l$range[1], l$range[2])
      if (l$invert) a <- 1 - a
      bg[[k]] <- .apply_opacity(f(a, bg[[k]]), bg[[k]], l$opacity)
    }
  }

  # NoData anywhere in the stack is NoData out - a blend of a hole is a hole
  for (k in seq_len(nb)) bg[[k]][is.na(tile)] <- NA
  bg
}

#' Render a blend stack to a raster
#'
#' Writes the composed image built by [rvt_blend()] as a Cloud-Optimized
#' GeoTIFF, values 0-1, ready to display - or straight to a WebP or JPEG
#' picture when `out_path` ends in `.webp` or `.jpg`.
#'
#' This is where the work happens. Layers are read tile by tile and blended in
#' one pass, so a stack of any depth costs a single traverse of the data.
#' Each layer's stretch range is measured once, before the first tile, which is
#' what keeps the result independent of `tile_size`.
#'
#' @section Output:
#' The file extension decides the format. `.tif` writes a Cloud-Optimized
#' GeoTIFF to use in a GIS. `.webp` and `.jpg` write a plain picture with no
#' georeferencing, for web pages and reports, with no GeoTIFF made along the
#' way; see [rvt_image()] for how the two formats compare. A single-band blend
#' is written in grey; for a colour palette use
#' `rvt_image(stack, "out.webp", col = ...)`.
#'
#' @param stack an `rvt_stack` from [rvt_blend()]
#' @param out_path where to write; defaults to a temporary GeoTIFF
#' @param quality for `.webp` and `.jpg`, compression quality 1-100 (default 90)
#' @inheritParams rvt_svf
#' @return `out_path`, invisibly
#' @seealso [rvt_blend()]
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_hillshade(dem) |>
#'   rvt_blend(rvt_svf(dem), "multiply", opacity = 0.25) |>
#'   rvt_render()
#' @export
rvt_render <- function(stack, out_path = fs::file_temp(ext = "tif"),
                        quality = 90, tile_size = NULL, threads = rvt_threads(),
                        overwrite = FALSE, progress = FALSE) {
  if (!inherits(stack, "rvt_stack"))
    stop("`stack` must come from rvt_blend() or rvt_stack()", call. = FALSE)
  .check_quality(quality)
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  format <- .image_format(out_path)

  r <- .render_setup(stack)
  on.exit(.rm_path(r$temp), add = TRUE)
  g <- stack$grid
  if (is.null(tile_size)) tile_size <- .auto_tile_size(g$nx, g$ny, 0L)

  .process_tiled(r$src, list(blend = out_path), 0L, tile_size,
                  function(tile, xres, yres, ctx) {
                    b <- .blend_tile(tile, r$layers, r$nb, ctx)
                    list(blend = if (is.null(format)) b else .image_bands(b, format))
                  },
                  progress = progress, threads = threads,
                  nbands = c(blend = if (is.null(format)) r$nb
                                     else .image_nbands(r$nb, format)),
                  aux = r$aux, image = !is.null(format), quality = quality)
  invisible(out_path)
}

## Everything a render needs before the first tile: ranges frozen, coarser
## layers upsampled, and every band of every layer as an aux entry. Shared by
## rvt_render() and rvt_image(). `temp` lists files the caller must remove.
.render_setup <- function(stack) {
  layers <- stack$layers
  # Freeze every range *now*, before a single tile is read. A range derived
  # per tile would make the same ground read differently depending on which
  # tile it fell in - the rvt-py sky-illumination bug this package refuses to
  # reproduce (see ?rvt_sky_illumination).
  layers <- lapply(layers, function(l) {
    if (is.null(l$range)) {
      # A multi-band layer gets ONE range spanning every band, not a range per
      # band. Stretching each channel to its own extremes would shift the
      # colour balance of an RGB image, and would rescale the three scales of
      # an rvt_mstp() result relative to each other. Stretch per channel
      # beforehand if that is what you actually want.
      rs <- vapply(seq_len(l$nbands),
                   function(b) rvt_range(l$path, band = b), numeric(2))
      l$range <- c(min(rs[1, ]), max(rs[2, ]))
    }
    if (!all(is.finite(l$range)) || l$range[1] >= l$range[2]) l$range <- c(0, 1)
    l
  })

  # Layers coarser than the finest are upsampled onto its grid - after the
  # ranges are frozen above, which read the original rasters.
  ups <- character(0)
  layers <- lapply(layers, function(l) {
    if (!.on_grid(l$path, stack$grid)) {
      l$path <- .upsample(l$path, stack$grid)
      ups <<- c(ups, l$path)
    }
    l
  })

  nb <- max(vapply(layers, function(l) l$nbands, integer(1)))
  if (nb > 1L && !all(vapply(layers, function(l) l$nbands, integer(1)) %in% c(1L, nb)))
    stop("layers must have 1 band or ", nb, " bands, not a mix of others",
          call. = FALSE)

  # every band of every layer becomes an aux entry; fac 1 and overlap 0 make
  # .process_tiled()'s window arithmetic reduce to the plain tile window
  aux <- list(); j <- 0L
  for (i in seq_along(layers)) {
    idx <- integer(0)
    for (b in seq_len(layers[[i]]$nbands)) {
      j <- j + 1L
      aux[[j]] <- list(path = layers[[i]]$path, band = b, fac = 1L, overlap = 0L)
      idx <- c(idx, j)
    }
    layers[[i]]$aux <- idx
  }

  list(layers = layers, aux = aux, nb = nb, src = layers[[1]]$path, temp = ups)
}
