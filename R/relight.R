## Relighting an orthophoto with the terrain beneath it.
##
## Not a blend mode, deliberately. The photo is multiplied by a *coloured
## light* - warm sun where the ground is lit, cool sky light everywhere, less
## of it in hollows - and that light runs above 1 on slopes facing the sun.
## rvt_blend() stretches every layer to 0-1 before blending, so it cannot
## express this; one function that computes the light and applies it can.

.relight_args <- c("sun_azimuth", "sun_elevation", "softness", "reach",
                   "saturation", "contrast", "vegetation", "sun_strength",
                   "sky_strength", "bounce", "occlusion", "sun_color",
                   "sky_color")

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
#'   \item{`verdant`}{Woodland and fields brought forward: deeper, more varied
#'     greens, with everything else left as it is.}
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
                 reach = 150, saturation = 1.12, contrast = 1.0,
                 vegetation = 0, sun_strength = 0.58, sky_strength = 0.66,
                 bounce = 0.74, occlusion = 0.20,
                 sun_color = "#FFEDCC", sky_color = "#A7BCE8"),
  neutral = list(sun_color = "#EEEEEE", sky_color = "#A5A5A5"),
  # the open-sky settings come along, or this ends up darker than
  # soft_shadows, which would be the wrong way round
  no_cast_shadows = list(bounce = 1, occlusion = 0.10, sky_strength = 0.70),
  soft_shadows = list(bounce = 0.86, sky_strength = 0.70, occlusion = 0.10),
  vivid = list(saturation = 1.26, contrast = 1.06),
  low_sun = list(sun_elevation = 15),
  high_sun = list(sun_elevation = 40),
  hard_shadows = list(softness = 0),
  verdant = list(vegetation = 0.55, saturation = 1.02,
                 sun_color = "#FFF3D6", sky_color = "#B0C6EA"),
  golden_hour = list(sun_color = "#FFE0B8", sky_color = "#8599D6",
                     saturation = 1.06),
  winter = list(sun_elevation = 18, sun_strength = 0.56,
                sun_color = "#EEF3FF", sky_color = "#A8C2F5",
                saturation = 0.90),
  # a 20 m rise at 6-8 degrees casts a shadow of 140-190 m, hence the reach
  twilight = list(sun_elevation = 8, softness = 3, reach = 300,
                  sun_strength = 0.45, sky_strength = 0.64, bounce = 0.60,
                  sun_color = "#F2B9A0", sky_color = "#8A86C4",
                  saturation = 1.0, contrast = 1.0),
  blue_hour = list(sun_elevation = 6, softness = 4, reach = 300,
                   sun_strength = 0.15, sky_strength = 0.88, bounce = 0.70,
                   occlusion = 0.50, sun_color = "#C9C2E6",
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
  if (p$vegetation < 0) stop("`vegetation` must not be negative.", call. = FALSE)
  p$sun_rgb <- as.numeric(grDevices::col2rgb(p$sun_color)) / 255
  p$sky_rgb <- as.numeric(grDevices::col2rgb(p$sky_color)) / 255
  p
}

## Colour of the photo before the light is applied: overall saturation, plus
## the optional vegetation treatment. Every step is per cell against fixed
## constants - never a tile statistic - so tiles stay independent.
.relight_colour <- function(photo, lum, p) {
  if (p$vegetation <= 0)
    return(lapply(photo, function(x) lum + (x - lum) * p$saturation))

  # how far the green channel leads the other two; 0.06 of full scale is
  # already unmistakably green, so foliage reaches full weight and grey roofs,
  # water and bare soil stay at 0
  w <- pmin(pmax((photo[[2]] - (photo[[1]] + photo[[3]]) / 2) / 0.06, 0), 1)
  ch <- lapply(photo, function(x) lum + (x - lum) * (p$saturation + p$vegetation * w))

  # a small fixed warm bias, not a stretch of the red-blue difference: the sky
  # light is blue, shaded canopy is already blue-leaning, and amplifying that
  # difference drove whole woods to teal instead of telling species apart
  shift <- 0.035 * p$vegetation * w
  ch[[1]] <- ch[[1]] + shift
  ch[[3]] <- ch[[3]] - shift

  # separation between light and dark foliage about a canopy-typical mid-tone,
  # which is where the nuance between stands comes from
  lapply(ch, function(x) 0.35 + (x - 0.35) * (1 + 0.30 * p$vegetation * w))
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

  # sky light, reduced in hollows that see less of the sky. The knee starts
  # well below open ground's sky-view factor: ordinary canopy sits near 0.4,
  # and cutting it off there left woodland almost black.
  open <- pmin(pmax((svf - 0.35) / 0.45, 0), 1)
  open[is.na(open)] <- 1
  ambient <- 1 - p$occlusion * (1 - open)

  lum <- 0.2126 * photo[[1]] + 0.7152 * photo[[2]] + 0.0722 * photo[[3]]
  ch <- .relight_colour(photo, lum, p)
  lapply(1:3, function(k) {
    light <- p$sun_strength * direct * p$sun_rgb[k] +
      p$sky_strength * ambient * p$sky_rgb[k]
    v <- 0.5 + (ch[[k]] * light - 0.5) * p$contrast
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
#' Give `out_path` a `.webp` or `.jpg` extension to write a plain picture
#' instead, for web pages and reports: no georeferencing, and no GeoTIFF made
#' along the way. WebP is the smaller of the two; JPEG is readable everywhere
#' and allows larger images. See [rvt_image()] for the comparison.
#'
#' @param ortho path to a three-band orthophoto
#' @param dem surface or terrain model covering the orthophoto's extent, or a
#'   vector of paths to mosaic. Its resolution may differ from the photo's by a
#'   whole factor: a 1 m surface under a 0.2 m photo is lit at 1 m and the
#'   light upsampled, so the result keeps the photo's full resolution.
#' @param out_path where to write; defaults to a temporary file
#' @param style name of a style in [rvt_relight_styles] (default `"default"`)
#' @param sun_azimuth,sun_elevation sun position in degrees
#' @param softness spread of the sun positions averaged for soft shadows, in
#'   degrees; 0 for hard shadows
#' @param reach how far shadows are traced, in map units
#' @param saturation colour strength of the photo; 1 leaves it unchanged
#' @param contrast contrast of the result; 1 leaves it unchanged
#' @param vegetation extra colour given to green ground only, 0 for none.
#'   Deepens foliage, pulls warm greens apart from cool ones and separates
#'   light from dark canopy, leaving roofs, roads, water and bare soil alone
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
                         vegetation = NULL,
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
  dem_grid <- .grid(dem)
  g <- .grid_join(dem_grid, .grid(ortho), "ortho", base = "`dem`")

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
  made <- c(hs, svf, shadows, scratch)
  on.exit(.rm_path(made), add = TRUE)

  # Whichever side is coarser is upsampled onto the finer grid. A 1 m surface
  # under a 0.2 m photo keeps every pixel of the photo: the light is computed
  # at the surface's own resolution, where it costs 25 times less, and is
  # smooth at 0.2 m anyway.
  if (dem_grid$nx != g$nx || dem_grid$ny != g$ny) {
    hs <- .upsample(hs, g); svf <- .upsample(svf, g)
    shadows <- vapply(shadows, .upsample, character(1), g = g, USE.NAMES = FALSE)
    made <- c(made, hs, svf, shadows)
  }
  if (!.on_grid(ortho, g)) {
    ortho <- .upsample(ortho, g)
    made <- c(made, ortho)
  }

  aux <- c(lapply(1:3, function(b)
             list(path = ortho, band = b, fac = 1L, overlap = 0L)),
           list(list(path = svf, fac = 1L, overlap = 0L)),
           lapply(shadows, function(s) list(path = s, fac = 1L, overlap = 0L)))
  info <- .dem_info(hs)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, 0L)

  # a picture: written straight from the tiled pass, no GeoTIFF in between
  format <- .image_format(out_path)
  if (!is.null(format)) {
    .process_tiled(hs, list(relight = out_path), 0L, tile_size,
                   function(tile, xres, yres, ctx) {
                     b <- .relight_tile(tile, ctx, p, length(shadows))
                     list(relight = .image_bands(lapply(b, `/`, 255), format))
                   },
                   progress = progress, threads = threads,
                   nbands = c(relight = .image_nbands(3L, format)), aux = aux,
                   image = TRUE, quality = quality)
    return(invisible(out_path))
  }

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
