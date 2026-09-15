## Relighting an orthophoto with the terrain beneath it.
##
## Not a blend mode, deliberately. The photo is multiplied by a *coloured
## light* - warm sun where the ground is lit, cool sky light everywhere, less
## of it in hollows - and that light runs above 1 on slopes facing the sun.
## rvt_blend() stretches every layer to 0-1 before blending, so it cannot
## express this; one function that computes the light and applies it can.

.relight_args <- c("sun_azimuth", "sun_elevation", "softness", "reach",
                   "saturation", "contrast", "sun_strength", "sky_strength",
                   "bounce", "occlusion", "sun_color", "sky_color")

#' Lighting styles for rvt_relight()
#'
#' Named sets of lighting settings for [rvt_relight()]. `default` holds every
#' setting; each other style lists only what it changes from `default`. Any
#' setting passed to `rvt_relight()` directly overrides the style.
#'
#' \describe{
#'   \item{`default`}{Warm afternoon sun from the north-west with soft blue
#'     shadows and slightly lifted colour.}
#'   \item{`neutral`}{The same light without colour: grey sun and sky, to see
#'     what the tints add.}
#'   \item{`no_cast_shadows`}{Shading and occlusion only; nothing casts a
#'     shadow onto its surroundings.}
#'   \item{`soft_shadows`}{Lighter, more open shadows.}
#'   \item{`vivid`}{Stronger colour and contrast.}
#'   \item{`low_sun`, `high_sun`}{A lower sun for long shadows, or a higher
#'     one for short shadows.}
#'   \item{`hard_shadows`}{Crisp shadow edges instead of soft ones.}
#'   \item{`golden_hour`}{Warm, low-contrast evening light.}
#'   \item{`winter`}{Cold, frosty light with longer shadows and dimmer
#'     highlights.}
#'   \item{`twilight`}{A very low rose-peach sun under a lavender sky.}
#'   \item{`blue_hour`}{The sun barely present; the scene lit mostly by a cool
#'     sky.}
#' }
#'
#' @format A named list of lists.
#' @seealso [rvt_relight()]
#' @examples
#' names(rvt_relight_styles)
#' rvt_relight_styles$golden_hour
#' @export
rvt_relight_styles <- list(
  default = list(sun_azimuth = 315, sun_elevation = 25, softness = 2,
                 reach = 150, saturation = 1.15, contrast = 1.05,
                 sun_strength = 0.68, sky_strength = 0.55,
                 bounce = 0.30, occlusion = 0.60,
                 sun_color = "#FFEDCC", sky_color = "#8CA6E6"),
  neutral = list(sun_color = "#EEEEEE", sky_color = "#A5A5A5"),
  no_cast_shadows = list(bounce = 1),
  soft_shadows = list(bounce = 0.45, sky_strength = 0.62),
  vivid = list(saturation = 1.30, contrast = 1.12),
  low_sun = list(sun_elevation = 15),
  high_sun = list(sun_elevation = 40),
  hard_shadows = list(softness = 0),
  golden_hour = list(sun_color = "#FFE0B8", sky_color = "#8599D6",
                     saturation = 1.08),
  winter = list(sun_elevation = 18, sun_strength = 0.60,
                sun_color = "#EEF3FF", sky_color = "#A8C2F5",
                saturation = 0.90),
  # a 20 m rise at 6-8 degrees casts a shadow of 140-190 m, hence the reach
  twilight = list(sun_elevation = 8, softness = 3, reach = 300,
                  sun_strength = 0.45, sky_strength = 0.60, bounce = 0.45,
                  sun_color = "#F2B9A0", sky_color = "#8A86C4",
                  saturation = 1.0, contrast = 1.0),
  blue_hour = list(sun_elevation = 6, softness = 4, reach = 300,
                   sun_strength = 0.15, sky_strength = 0.85, bounce = 0.60,
                   occlusion = 0.70, sun_color = "#C9C2E6",
                   sky_color = "#6C86CC", saturation = 0.92, contrast = 1.0)
)

## The style's settings with any explicitly given ones laid over them.
.relight_params <- function(style, overrides) {
  if (!is.character(style) || length(style) != 1L ||
      !style %in% names(rvt_relight_styles))
    stop(sprintf("Unknown style. Choose one of: %s.",
                 paste(names(rvt_relight_styles), collapse = ", ")),
         call. = FALSE)
  p <- utils::modifyList(rvt_relight_styles$default, rvt_relight_styles[[style]])
  overrides <- overrides[!vapply(overrides, is.null, logical(1))]
  p <- utils::modifyList(p, overrides)

  if (p$sun_elevation - p$softness / 2 <= 0 || p$sun_elevation >= 90)
    stop("`sun_elevation` must stay between 0 and 90 degrees, softness included.",
         call. = FALSE)
  if (p$softness < 0) stop("`softness` must not be negative.", call. = FALSE)
  p$sun_rgb <- as.numeric(grDevices::col2rgb(p$sun_color)) / 255
  p$sky_rgb <- as.numeric(grDevices::col2rgb(p$sky_color)) / 255
  p
}

## One tile. `tile` is the hillshade; the aux windows arrive as the three
## orthophoto bands, the sky-view factor, then the shadow runs.
.relight_tile <- function(tile, ctx, p, n_shadow) {
  a <- ctx$aux
  photo <- lapply(1:3, function(k) a[[k]]$m / 255)
  svf <- a[[4]]$m
  shadow <- Reduce(`+`, lapply(a[4L + seq_len(n_shadow)], function(w) w$m)) / n_shadow

  # direct light: hillshade scaled so flat open ground is 1, dimmed where cast
  # shadow falls but never to zero - bounce light fills shadows in reality
  shade <- p$bounce + (1 - p$bounce) * shadow
  direct <- tile * shade / sin(p$sun_elevation * pi / 180)
  direct[is.na(direct)] <- 1
  direct <- pmin(direct, 1.6)

  # sky light, reduced in hollows that see less of the sky
  open <- pmin(pmax((svf - 0.6) / 0.4, 0), 1)
  open[is.na(open)] <- 1
  ambient <- 1 - p$occlusion * (1 - open)

  lum <- 0.2126 * photo[[1]] + 0.7152 * photo[[2]] + 0.0722 * photo[[3]]
  lapply(1:3, function(k) {
    ch <- lum + (photo[[k]] - lum) * p$saturation
    light <- p$sun_strength * direct * p$sun_rgb[k] +
      p$sky_strength * ambient * p$sky_rgb[k]
    v <- 0.5 + (ch * light - 0.5) * p$contrast
    pmin(pmax(v, 0), 1) * 255
  })
}

#' Relight an orthophoto with its terrain
#'
#' Lights an aerial photo with the terrain beneath it, the way a rendered
#' landscape model is lit: a warm sun on slopes that face it, soft cast
#' shadows, cool sky light that colours the shadows, and less of that sky
#' light in enclosed ground. The result reads as a small physical model rather
#' than a flat photo with shading laid over it.
#'
#' Use a **surface model** for `dem` when the photo shows trees and buildings,
#' so that they cast shadows too, and a terrain model when you want only the
#' shape of the ground.
#'
#' @section How it works:
#' Hillshade, sky-view factor and cast shadow are computed from `dem`. Shadows
#' are averaged over sun positions `softness` degrees apart, which softens
#' their edges. Light is then built per cell as sun light plus sky light:
#'
#' * **sun light**: `sun_color` times `sun_strength`, scaled by how squarely
#'   the ground faces the sun, and reduced in cast shadow to a fraction
#'   `bounce` of full strength rather than to nothing;
#' * **sky light**: `sky_color` times `sky_strength`, reduced by up to
#'   `occlusion` where the sky-view factor says the ground is enclosed.
#'
#' The photo's colours are strengthened by `saturation`, multiplied by that
#' light, and given `contrast`. Because the light is coloured, shadows pick up
#' a cool cast from the sky light; the tint is subtle in deep shadow, where
#' little light of any colour remains.
#'
#' Sunlit open ground keeps close to the photo's brightness, but anything in
#' shadow is darker, so heavily wooded or built-up areas come out noticeably
#' darker overall. Raise `bounce` or lower `occlusion` to lighten them.
#'
#' @section Styles:
#' `style` picks a named set of settings from [rvt_relight_styles], such as
#' `"golden_hour"` or `"winter"`. Any setting given here overrides the style,
#' so `style = "winter", sun_elevation = 10` is a winter light with an even
#' lower sun.
#'
#' @section Output:
#' A three-band, 8-bit Cloud-Optimized GeoTIFF compressed with WebP. `quality`
#' trades file size for fidelity; the default of 90 is visually lossless for
#' imagery at typical viewing sizes.
#'
#' @param ortho path to a three-band orthophoto
#' @param dem surface or terrain model on exactly the orthophoto's grid, or a
#'   vector of paths to mosaic
#' @param out_path where to write; defaults to a temporary file
#' @param style name of a style in [rvt_relight_styles] (default `"default"`)
#' @param sun_azimuth,sun_elevation sun position in degrees
#' @param softness spread of the sun positions averaged for soft shadows, in
#'   degrees; 0 for hard shadows
#' @param reach how far shadows are traced, in map units
#' @param saturation colour strength of the photo; 1 leaves it unchanged
#' @param contrast contrast of the result; 1 leaves it unchanged
#' @param sun_strength,sky_strength strength of the sun and sky light
#' @param bounce share of sun light that still reaches ground in cast shadow
#' @param occlusion how much sky light enclosed ground loses, 0 to 1
#' @param sun_color,sky_color colours of the sun and sky light
#' @param quality WebP quality, 1 to 100 (default 90)
#' @inheritParams rvt_svf
#' @return `out_path`, invisibly
#' @seealso [rvt_relight_styles], [rvt_blend()] for general layer blending,
#'   [rvt_data_lgln()] for sample data
#' @examples
#' \donttest{
#' pt <- c(9.9464, 51.6317)                      # Burg Hardenberg
#' rgb <- rvt_data_lgln(pt, "rgb", res = 1)
#' dsm <- rvt_data_lgln(pt, "dsm")
#' rvt_relight(rgb, dsm) |> rvt_plot(minmax_def = c(0, 0, 0, 255, 255, 255))
#' rvt_relight(rgb, dsm, style = "golden_hour") |>
#'   rvt_plot(minmax_def = c(0, 0, 0, 255, 255, 255))
#' }
#' @export
rvt_relight <- function(ortho, dem, out_path = fs::file_temp(ext = "tif"),
                         style = "default",
                         sun_azimuth = NULL, sun_elevation = NULL,
                         softness = NULL, reach = NULL,
                         saturation = NULL, contrast = NULL,
                         sun_strength = NULL, sky_strength = NULL,
                         bounce = NULL, occlusion = NULL,
                         sun_color = NULL, sky_color = NULL,
                         quality = 90, tile_size = NULL,
                         threads = rvt_threads(), overwrite = FALSE,
                         progress = FALSE) {
  p <- .relight_params(style, mget(.relight_args, envir = environment()))
  if (!is.numeric(quality) || length(quality) != 1L || quality < 1 || quality > 100)
    stop("`quality` must be a number from 1 to 100.", call. = FALSE)

  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  ortho <- .as_path(ortho)
  dem <- rvt_mosaic(dem)
  if (.nbands(ortho) != 3L)
    stop("`ortho` must have three bands (red, green, blue).", call. = FALSE)
  .check_aligned(dem, ortho, "ortho")

  tmp <- function() fs::file_temp(ext = "tif")
  hs <- rvt_hillshade(dem, tmp(), sun_azimuth = p$sun_azimuth,
                      sun_elevation = p$sun_elevation, threads = threads)
  svf <- rvt_svf(dem, tmp(), reach = 10, threads = threads)
  offsets <- if (p$softness == 0) 0 else c(-p$softness, 0, p$softness)
  shadows <- vapply(offsets, function(d)
    as.character(rvt_shadow(dem, tmp(),
                            sun_azimuth = (p$sun_azimuth + d) %% 360,
                            sun_elevation = p$sun_elevation + d / 2,
                            reach = p$reach, threads = threads)),
    character(1))
  scratch <- tmp()
  on.exit(.rm_path(c(hs, svf, shadows, scratch)), add = TRUE)

  aux <- c(lapply(1:3, function(b)
             list(path = ortho, band = b, fac = 1L, overlap = 0L)),
           list(list(path = svf, fac = 1L, overlap = 0L)),
           lapply(shadows, function(s) list(path = s, fac = 1L, overlap = 0L)))
  info <- .dem_info(hs)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, 0L)

  .process_tiled(hs, list(relight = scratch), 0L, tile_size,
                 function(tile, xres, yres, ctx)
                   list(relight = .relight_tile(tile, ctx, p, length(shadows))),
                 progress = progress, threads = threads,
                 nbands = c(relight = 3L), aux = aux)

  # 8-bit WebP inside a COG. No predictor: WebP does its own prediction, and
  # predictors only help general-purpose codecs such as DEFLATE.
  gdalraster::translate(scratch, out_path,
    cl_arg = c("-ot", "Byte", "-a_nodata", "none", "-of", "COG",
               "-co", "COMPRESS=WEBP",
               "-co", paste0("QUALITY=", round(quality)),
               "-co", "RESAMPLING=CUBIC",
               "-co", paste0("NUM_THREADS=", threads)),
    quiet = TRUE)
  invisible(out_path)
}
