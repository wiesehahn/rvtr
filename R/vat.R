## Visualization for Archaeological Topography (Kokalj & Somrak 2019).
##
## VAT stacks four normalised layers - hillshade at the bottom, then slope,
## positive openness and sky-view factor blended over it - and the "combined"
## product averages two such stacks, one tuned for general terrain and one for
## flat terrain. Layer parameters and blend modes below follow RVT's
## blender_VAT.json and default_terrains_settings.json; the stacking order and
## blend arithmetic follow rvt.blend.BlenderCombination$render_all_images(),
## which renders from the last layer to the first.

## RVT's normalize_lin(): linear cut-off then stretch to 0-1.
.norm_lin <- function(x, lo, hi) {
  x[x < lo] <- lo
  x[x > hi] <- hi
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

#' VAT terrain presets
#'
#' Parameter sets for [rvt_vat()], from RVT's `default_terrains_settings.json`.
#' `rvt_preset_general` suits typical relief; `rvt_preset_flat` uses a larger
#' search radius and more noise removal for subtle micro-relief on flat ground.
#' Pass a customised copy (e.g. `modifyList(rvt_preset_general, list(radius_max = 15))`)
#' to tune either preset.
#'
#' @format A list with `sun_elevation`, `slope`, `svf`, `opns` (normalisation
#'   ranges), `radius_max` and `noise_removal`.
#' @name rvt_presets
NULL

#' @rdname rvt_presets
#' @export
rvt_preset_general <- list(sun_elevation = 35, slope = c(0, 50), svf = c(0.7, 1),
                            opns = c(68, 93), radius_max = 10, noise_removal = 0)

#' @rdname rvt_presets
#' @export
rvt_preset_flat <- list(sun_elevation = 15, slope = c(0, 15), svf = c(0.9, 1),
                         opns = c(85, 93), radius_max = 20, noise_removal = 3)

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
  sa <- slope_hillshade(.pad_edge(tile, pad_sa), pad_sa, nrow(tile), ncol(tile),
                         xres, yres, 315,
                         vapply(presets, function(p) p$sun_elevation, numeric(1)),
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
                           TRUE, TRUE, threads)
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
#' Visualization for Archaeological Topography: a blend of hillshade, slope,
#' positive openness and sky-view factor, rendered once with a general-terrain
#' preset and once with a flat-terrain preset and averaged. Values are 0-1.
#'
#' @inheritParams rvt_svf
#' @param presets list of two parameter sets, general first then flat; see
#'   [rvt_presets]
#' @param rvt_compat reproduce a bug in rvt-py's blend_overlay() that makes
#'   the openness layer's declared 50% opacity have no effect. `FALSE`
#'   (default) applies the opacity as configured; `TRUE` is only useful for
#'   comparing output against rvt-py.
#' @return `out_path`, invisibly. `dem` is the first argument, so this is
#'   pipe-friendly: `dem |> rvt_vat("vat.tif")`.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_vat(dem)
#' @export
rvt_vat <- function(dem, out_path = tempfile(fileext = ".tif"),
                num_directions = 16,
                presets = list(rvt_preset_general, rvt_preset_flat),
                tile_size = NULL, threads = rvt_threads(),
                rvt_compat = FALSE, overwrite = FALSE, progress = FALSE) {
  if (!overwrite && file.exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  overlap <- as.integer(max(vapply(presets, function(p) p$radius_max, numeric(1))) + 1L)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

  .process_tiled(dem, list(vat = out_path), overlap, tile_size,
                  function(tile, xres, yres)
                    .vat_tile(tile, xres, yres, num_directions, presets, threads,
                               rvt_compat),
                  progress = progress, threads = threads)
  invisible(out_path)
}
