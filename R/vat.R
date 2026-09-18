## Visualization for Archaeological Topography (Kokalj & Somrak 2019).
##
## VAT stacks four normalised layers - hillshade at the bottom, then slope,
## positive openness and sky-view factor blended over it - and the "combined"
## product averages two such stacks, one tuned for general terrain and one for
## flat terrain. Layer parameters and blend modes below follow RVT's
## blender_VAT.json and default_terrains_settings.json; the stacking order and
## blend arithmetic follow rvt.blend.BlenderCombination$render_all_images(),
## which renders from the last layer to the first.

## `.norm_lin()`, `.blend_overlay()` and `.apply_opacity()` now live in
## R/blend.R, which generalises this stack into rvt_blend(). They behave
## identically; rvt_vat() is deliberately *not* reimplemented on top of the
## engine, because it is validated against rvt-py down to the rvt_compat bug
## flag and that is worth more than sharing the code path. A test asserts the
## engine reproduces .vat_stack() instead.

#' VAT terrain presets
#'
#' The two parameter sets that [rvt_vat()] renders and averages.
#' `rvt_preset_general` suits ordinary relief; `rvt_preset_flat` looks further
#' out and filters noise harder, to bring up the subtle micro-relief typical of
#' flat ground.
#'
#' Both are plain lists, so the way to adjust one is a modified copy:
#' `modifyList(rvt_preset_general, list(reach = 20))`.
#'
#' @format A list with:
#' \describe{
#'   \item{`sun_elevation`}{height of the sun above the horizon, in degrees,
#'     for the hillshade layer. Lower angles throw longer shadows and
#'     exaggerate faint relief - hence 15 for flat terrain against 35 for
#'     general.}
#'   \item{`slope`}{the slope range, in degrees, stretched across the layer's
#'     full brightness range. Slopes past the upper end are clipped, so a
#'     narrower range gives more contrast among gentle slopes (15 for flat
#'     terrain, 50 for general).}
#'   \item{`svf`}{the same idea for the sky-view factor layer - the
#'     sky-fraction range that gets stretched. `c(0.9, 1)` on flat terrain
#'     spreads out what would otherwise be an almost invisible band of
#'     values.}
#'   \item{`opns`}{the same again for the positive openness layer, in
#'     degrees.}
#'   \item{`reach`}{how far the horizon search behind the openness and
#'     sky-view factor layers looks, in map units (metres, normally).}
#'   \item{`noise_removal`}{0-3, ignoring progressively more of the innermost
#'     part of each search ray.}
#' }
#' @seealso [rvt_vat()]
#' @name rvt_presets
NULL

#' @rdname rvt_presets
#' @export
rvt_preset_general <- list(sun_elevation = 35, slope = c(0, 50), svf = c(0.7, 1),
                            opns = c(68, 93), reach = 10, noise_removal = 0)

#' @rdname rvt_presets
#' @export
rvt_preset_flat <- list(sun_elevation = 15, slope = c(0, 15), svf = c(0.9, 1),
                         opns = c(85, 93), reach = 20, noise_removal = 3)

## One VAT stack from already-computed ingredients.
##
## The openness layer blends over the hillshade+slope background in Overlay
## mode at a declared 50% opacity. rvt-py has a bug where its blend_overlay()
## mutates its `background` argument in place and returns it, so the opacity
## step that should follow ends up blending the result with itself and the
## declared 50% has no effect - the openness layer is applied at full
## strength. `rvt_compat = TRUE` reproduces that bug, which is only useful for
## checking output against rvt-py; the default applies the 50% opacity as
## configured.
.vat_stack <- function(hillshade, slope_deg, opns, svf_arr, preset, rvt_compat = FALSE) {
  # layer 4 (bottom): hillshade, Normal, opacity 100
  rendered <- .norm_lin(hillshade, 0, 1)

  # layer 3: slope, Luminosity, opacity 50. RVT inverts the slope scale, and
  # for single-band images blend_luminosity() reduces to the active layer.
  active <- 1 - .norm_lin(slope_deg, preset$slope[1], preset$slope[2])
  rendered <- .apply_opacity(active, rendered, 0.5)

  # layer 2: positive openness, Overlay, opacity 50 (see note above)
  active <- .norm_lin(opns, preset$opns[1], preset$opns[2])
  top <- .blend_overlay(active, rendered)
  rendered <- if (rvt_compat) top else .apply_opacity(top, rendered, 0.5)

  # layer 1 (top): sky-view factor, Multiply, opacity 25
  active <- .norm_lin(svf_arr, preset$svf[1], preset$svf[2])
  rendered <- .apply_opacity(active * rendered, rendered, 0.25)

  rendered
}

.vat_tile <- function(tile, xres, yres, num_directions, presets, threads,
                       rvt_compat = FALSE) {
  # edge padding for the derivatives, matching rvt.vis.slope_aspect()
  pad_sa <- 1L
  sun_el <- vapply(presets, function(p) p$sun_elevation, numeric(1))
  sa <- slope_hillshade(.pad_edge(tile, pad_sa), pad_sa, nrow(tile), ncol(tile),
                         xres, yres, rep(315, length(sun_el)), sun_el,
                         threads)
  slope_deg <- sa$slope * 180 / pi

  stacks <- vector("list", length(presets))
  for (i in seq_along(presets)) {
    p <- presets[[i]]
    off <- .direction_offsets(num_directions, p$radius_max,
                               .radius_min(p$radius_max, p$noise_removal), xres, yres)
    pad <- as.integer(p$radius_max + 1L)
    h <- horizon_svf_opns(.pad_reflect(tile, pad), pad, nrow(tile), ncol(tile),
                           off$dx, off$dy, off$dist, off$starts, off$ends,
                           TRUE, TRUE, FALSE, numeric(0), threads)
    stacks[[i]] <- .vat_stack(sa$hillshade[[i]], slope_deg, h$opns, h$svf, p, rvt_compat)
  }

  # "combined": the general stack laid over the flat stack at 50% opacity,
  # i.e. a straight average of the two.
  combined <- .apply_opacity(stacks[[1]], stacks[[2]], 0.5)
  combined[is.na(tile)] <- NA
  list(vat = combined)
}

#' Archaeological VAT (combined)
#'
#' Visualization for Archaeological Topography: a ready-made composite for
#' spotting earthworks, combining four different views of the terrain into a
#' single greyscale image so you don't have to read several rasters side by
#' side. Values are 0-1, ready to display directly.
#'
#' Each ingredient contributes something the others lack: hillshade gives the
#' eye the familiar sense of relief, slope sharpens edges, positive openness
#' lifts convex features such as banks and mounds, and sky-view factor darkens
#' enclosed ground such as ditches. The whole thing is rendered twice - once
#' with settings for ordinary relief, once tuned for flat ground where
#' features are subtle - and the two are averaged, so a single run works
#' reasonably across both without retuning.
#'
#' If you want to interpret actual numbers rather than look at a picture, this
#' is the wrong tool: it is a display product, and the blending discards the
#' physical units of its ingredients. Use [rvt_msrm()] (heights in metres) or
#' [rvt_openness()] (angles in degrees) for that.
#'
#' @section How it works:
#' For each of the two presets ([rvt_presets]), the DEM is turned into four
#' layers, each stretched to 0-1 over a range set by the preset, then stacked
#' from the bottom up:
#'
#' * **hillshade** at the preset's sun elevation, from the north-west - the
#'   background everything else is laid over;
#' * **slope**, inverted so flat ground is bright, blended over it at half
#'   strength;
#' * **positive openness**, blended in Overlay mode, which lightens what is
#'   already light and darkens what is already dark, so convex features gain
#'   contrast without washing the image out;
#' * **sky-view factor**, multiplied in at a quarter strength, which only
#'   darkens - deepening enclosed ground.
#'
#' The two resulting images are then averaged. The horizon searches use each
#' preset's own `reach` and `noise_removal`, so the flat-terrain pass
#' looks further out and filters more aggressively than the general one.
#' Computation is at native resolution, tile by tile.
#'
#' @section Tuning:
#' * The two presets are ordinary lists, so the easiest adjustment is a
#'   modified copy - for example
#'   `modifyList(rvt_preset_general, list(reach = 20))` to make the
#'   general pass respond to larger features. See [rvt_presets] for the
#'   fields and what they mean.
#' * To render a single terrain type rather than the average of two, pass the
#'   same preset twice, e.g. `presets = list(rvt_preset_flat, rvt_preset_flat)`.
#'   This is worth doing when your area really is uniformly flat, since the
#'   general-terrain pass otherwise dilutes the contrast the flat preset was
#'   tuned to give.
#' * `num_directions` is shared by both passes and behaves as in [rvt_svf()].
#' * `rvt_compat` exists only for checking output against rvt-py and should
#'   be left alone otherwise; see its parameter description.
#'
#' @inheritParams rvt_svf
#' @param presets list of two parameter sets, general first then flat; see
#'   [rvt_presets]
#' @param rvt_compat reproduce a bug in rvt-py's blend_overlay() that makes
#'   the openness layer's declared 50% opacity have no effect. `FALSE`
#'   (default) applies the opacity as configured; `TRUE` is only useful for
#'   comparing output against rvt-py.
#' @param quality for a `.webp` or `.jpg` `out_path`, compression quality 1-100
#'   (default 90)
#' @return `out_path`, invisibly. `dem` is the first argument, so this is
#'   pipe-friendly: `dem |> rvt_vat("vat.tif")`.
#' @section Output:
#' A Cloud-Optimized GeoTIFF, values 0-1. Give `out_path` a `.webp` or `.jpg`
#' extension to write a plain greyscale picture instead, with no GeoTIFF made
#' along the way and no georeferencing; see [rvt_image()] for the formats.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' dem |> rvt_vat()
#' @export
rvt_vat <- function(dem, out_path = fs::file_temp(ext = "tif"),
                num_directions = 16,
                presets = list(rvt_preset_general, rvt_preset_flat),
                tile_size = NULL, threads = rvt_threads(),
                rvt_compat = FALSE, quality = 90, overwrite = FALSE,
                progress = FALSE) {
  .check_quality(quality)
  out_path <- .out_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)
  format <- .image_format(out_path)

  info <- .dem_info(dem)
  # presets carry `reach` in map units; the kernel wants whole cells
  presets <- lapply(presets, function(p) {
    p$radius_max <- .cells(p$reach, info$xres, "reach", "outer")
    p
  })
  overlap <- as.integer(max(vapply(presets, function(p) p$radius_max, numeric(1))) + 1L)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

  .process_tiled(dem, list(vat = out_path), overlap, tile_size,
                  function(tile, xres, yres) {
                    v <- .vat_tile(tile, xres, yres, num_directions, presets,
                                   threads, rvt_compat)
                    if (!is.null(format)) v$vat <- .image_bands(list(v$vat), format)
                    v
                  },
                  progress = progress, threads = threads,
                  nbands = if (!is.null(format)) c(vat = .image_nbands(1L, format)),
                  image = !is.null(format), quality = quality)
  invisible(out_path)
}
