## Display images: WebP and JPEG, written straight from the tiled pass.
##
## GDAL's WEBP and JPEG drivers are CreateCopy-only - they cannot be written
## tile by tile - so the tiled pass fills its usual scratch GeoTIFF and a single
## translate turns that into the image. No COG is made along the way.
##
## These are pictures for looking at, not data: visually lossless is the bar,
## and every format here is lossy at a fixed quality. Anything that needs exact
## values stays a GeoTIFF.
##
## Two formats, chosen by measurement against PNG and lossless WebP on a relit
## 0.2 m orthophoto, a hillshade blend and a geomorphon class map:
## * WebP is the smallest (photo 5.8 MB at q90 against 6.3 MB JPEG, 41 MB PNG,
##   60 MB DEFLATE COG) and has an alpha channel for NoData.
## * JPEG for compatibility, images up to 65,500 px (WebP stops at 16,383) and
##   speed: 0.8 s against 8.1 s for the 5000 x 5000 photo.
## PNG was dropped as far larger. Lossless WebP was considered for class maps,
## where lossy compression shifts the colour of single-cell classes, and
## decided against: at viewing size that is not visible, and it is not what
## these images are for.

.image_drivers <- c(webp = "WEBP", jpg = "JPEG", jpeg = "JPEG")
.image_limits <- c(WEBP = 16383, JPEG = 65500)

## The image driver an output path asks for, or NULL for a GeoTIFF.
.image_format <- function(path) {
  ext <- tolower(fs::path_ext(path))
  if (ext %in% names(.image_drivers)) unname(.image_drivers[ext]) else NULL
}

## Bands written for `n` input bands (1 or 3): WebP always gets RGB plus an
## alpha band, so NoData is transparent; JPEG keeps 1 or 3 and has no alpha.
.image_nbands <- function(n, format) if (format == "WEBP") 4L else as.integer(n)

## A list of 0-1 matrices (1 or 3) -> 0-255 matrices ready to write.
.image_bands <- function(bands01, format) {
  hole <- Reduce(`|`, lapply(bands01, is.na))
  v <- lapply(bands01, function(m) {
    m <- round(pmin(pmax(m, 0), 1) * 255)
    m[hole] <- 0
    m
  })
  if (format == "WEBP") {
    if (length(v) == 1L) v <- rep(v, 3L)
    alpha <- matrix(255, nrow(hole), ncol(hole))
    alpha[hole] <- 0
    v <- c(v, list(alpha))
  }
  v
}

.image_check_size <- function(nx, ny, format) {
  lim <- .image_limits[[format]]
  if (nx <= lim && ny <= lim) return(invisible(TRUE))
  alt <- if (format == "WEBP") " JPEG (`.jpg`) allows up to 65,500 px;" else ""
  stop(sprintf(paste0(
    "This image would be %d x %d px, but %s allows at most %s px per side.%s ",
    "use a coarser resolution (`rvt_resample()`), a smaller area, or `.tif`."),
    nx, ny, if (format == "WEBP") "WebP" else "JPEG",
    format(lim, big.mark = ","), alt), call. = FALSE)
}

.check_quality <- function(quality) {
  if (!is.numeric(quality) || length(quality) != 1L || is.na(quality) ||
      quality < 1 || quality > 100)
    stop("`quality` must be a number from 1 to 100.", call. = FALSE)
  invisible(TRUE)
}

## Scratch GeoTIFF holding 0-255 -> the image. PAM is switched off so GDAL
## writes no .aux.xml beside it: the images are plain pictures, by decision.
.finalize_image <- function(scratch_path, out_path, format, quality) {
  old <- gdalraster::get_config_option("GDAL_PAM_ENABLED")
  gdalraster::set_config_option("GDAL_PAM_ENABLED", "NO")
  on.exit(gdalraster::set_config_option("GDAL_PAM_ENABLED", old), add = TRUE)
  .rm_path(out_path)
  gdalraster::translate(scratch_path, out_path, quiet = TRUE,
    cl_arg = c("-of", format, "-ot", "Byte", "-a_nodata", "none",
               "-co", paste0("QUALITY=", round(quality))))
  .rm_path(scratch_path)
  invisible(out_path)
}

## 256 colours as a 3 x 256 matrix, from anything rvt_plot()'s `col` accepts.
.col_table <- function(col) {
  cols <- .resolve_col(col, 256L)
  cols <- if (is.function(cols)) cols(256L)
          else if (length(cols) != 256L) grDevices::colorRampPalette(cols)(256L)
          else cols
  grDevices::col2rgb(cols) / 255
}

## One band stretched to 0-1 -> three 0-1 bands through a colour table.
.apply_col <- function(v01, table) {
  idx <- pmin(floor(v01 * 256), 255) + 1
  lapply(1:3, function(k) {
    out <- table[k, idx]
    dim(out) <- dim(v01)
    out
  })
}

#' Save a raster or blend as an image
#'
#' Writes a WebP or JPEG picture of any raster or blend stack, at full
#' resolution: a colour palette and stretch for a single band, the colours as
#' they are for three bands. Where [rvt_plot()] draws a quick look on screen,
#' this makes the file you put on a web page or in a report.
#'
#' The format follows the file extension: `.webp` or `.jpg`. Images carry no
#' georeferencing; for a result to use in a GIS, keep the `.tif` the metric
#' function wrote.
#'
#' [rvt_render()], [rvt_relight()] and [rvt_vat()] write images directly too -
#' give them an image path - so their results need no second step.
#'
#' @section Choosing a format:
#' Both are compressed for viewing, not for keeping values: the default quality
#' of 90 looks the same as the original, but the numbers are not preserved.
#' Keep the GeoTIFF for anything that needs them.
#'
#' **WebP** is the smaller file in every case measured - about 10% smaller than
#' JPEG for photographs and shading, and a tenth of a compressed GeoTIFF - and
#' shows NoData as transparent.
#'
#' **JPEG** is the one to reach for when a tool cannot read WebP, or for images
#' wider than WebP's limit of 16,383 pixels (JPEG allows 65,500). It has no
#' transparency, so NoData is black.
#'
#' @section Stretch and colours:
#' A single band is stretched across `range` and coloured with `col`, exactly
#' as in [rvt_plot()]. `range = NULL` uses the raster's own minimum and
#' maximum; [rvt_range()] gives a percentile stretch, and each metric's help
#' page lists recommended ranges. The range is fixed before any tile is read,
#' so the picture does not depend on `tile_size`.
#'
#' Three bands - an RGB result such as [rvt_mstp()], or a blend stack over an
#' orthophoto - are stretched over one range spanning all bands, and `col` is
#' not used. A blend stack is stretched layer by layer when it is blended (see
#' [rvt_blend()]), so `range` does not apply to it either.
#'
#' @param x a raster path, or an `rvt_stack` from [rvt_blend()]
#' @param out_path the image to write, ending in `.webp` or `.jpg`
#' @param col palette name or colours for a single band, as in [rvt_plot()]
#'   (default `"Grays"`); see [rvt_palettes()]
#' @param range `c(lo, hi)` to stretch a single-band or RGB raster over, or
#'   `NULL` (default) for its own minimum and maximum
#' @param band which band to draw, or three band numbers for red, green and
#'   blue. Defaults to band 1, or 1-3 for a three-band raster.
#' @param quality compression quality, 1-100 (default 90)
#' @inheritParams rvt_svf
#' @return `out_path`, invisibly
#' @seealso [rvt_plot()] for a quick look on screen, [rvt_render()] for blends
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' svf <- rvt_svf(dem)
#' rvt_image(svf, tempfile(fileext = ".webp"), range = c(0.7, 1))
#'
#' # signed data: a diverging palette with a symmetric range
#' dem |> rvt_slrm() |>
#'   rvt_image(tempfile(fileext = ".webp"), "Blue-Red 3", range = c(-1, 1))
#' @export
rvt_image <- function(x, out_path, col = "Grays", range = NULL, band = NULL,
                      quality = 90, tile_size = NULL, threads = rvt_threads(),
                      overwrite = FALSE, progress = FALSE) {
  out_path <- .as_path(out_path)
  format <- .image_format(out_path)
  if (is.null(format))
    stop("`out_path` must end in `.webp` or `.jpg`.", call. = FALSE)
  .check_quality(quality)
  if (!is.null(range) && (length(range) != 2L || !all(is.finite(range)) ||
                           range[1] >= range[2]))
    stop("`range` must be `c(lo, hi)` with lo < hi, or NULL", call. = FALSE)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  table <- .col_table(col)

  if (inherits(x, "rvt_stack")) {
    r <- .render_setup(x)
    on.exit(.rm_path(r$temp), add = TRUE)
    g <- x$grid
    if (is.null(tile_size)) tile_size <- .auto_tile_size(g$nx, g$ny, 0L)
    .process_tiled(r$src, list(image = out_path), 0L, tile_size,
      function(tile, xres, yres, ctx) {
        b <- .blend_tile(tile, r$layers, r$nb, ctx)
        if (r$nb == 1L) b <- .apply_col(b[[1]], table)
        list(image = .image_bands(b, format))
      },
      progress = progress, threads = threads, aux = r$aux,
      nbands = c(image = .image_nbands(3L, format)),
      image = TRUE, quality = quality)
    return(invisible(out_path))
  }

  x <- .as_path(x)
  nb_src <- .nbands(x)
  if (is.null(band)) band <- if (nb_src == 3L) 1:3 else 1L
  if (!length(band) %in% c(1L, 3L) || any(band < 1) || any(band > nb_src))
    stop(sprintf("`band` must be one or three band numbers from 1 to %d.", nb_src),
         call. = FALSE)
  if (is.null(range)) {
    rs <- vapply(band, function(b) rvt_range(x, band = b), numeric(2))
    range <- c(min(rs[1, ]), max(rs[2, ]))
    if (!all(is.finite(range)) || range[1] >= range[2]) range <- c(range[1], range[1] + 1)
  }

  info <- .dem_info(x, band[1])
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, 0L)
  aux <- if (length(band) == 3L)
    lapply(band[2:3], function(b) list(path = x, band = b, fac = 1L, overlap = 0L))

  .process_tiled(x, list(image = out_path), 0L, tile_size,
    function(tile, xres, yres, ctx = NULL) {
      s <- function(m) .norm_lin(m, range[1], range[2])
      b <- if (length(band) == 1L) .apply_col(s(tile), table)
           else c(list(s(tile)), lapply(ctx$aux, function(w) s(w$m)))
      list(image = .image_bands(b, format))
    },
    band = band[1], progress = progress, threads = threads, aux = aux,
    nbands = c(image = .image_nbands(3L, format)),
    image = TRUE, quality = quality)
  invisible(out_path)
}
