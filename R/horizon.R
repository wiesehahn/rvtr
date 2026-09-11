## Sky-view factor (isotropic and anisotropic), positive and negative openness.

## Per-direction weights for the anisotropic sky-view factor: 1 towards
## `main_direction`, falling to `min_weight` at the opposite azimuth.
##
## The directions produced by .direction_offsets() run counterclockwise from
## east (dx = cos, dy = -sin, and the row index grows southward), so direction
## d points at compass azimuth 90 - angle_d. Converting here rather than
## weighting the raw index is what makes `main_direction = 315` actually mean
## "brightest towards the north-west", matching the sun azimuth convention
## used for hillshade.
.asvf_weights <- function(num_directions, main_direction, poly_level, min_weight) {
  angle_deg <- (360 / num_directions) * (0:(num_directions - 1))
  azimuth <- (90 - angle_deg) %% 360
  delta <- (azimuth - main_direction) * pi / 180
  (1 - min_weight) * cos(delta / 2)^poly_level + min_weight
}

## Level 1 is mild anisotropy, level 2 strong (the opposite sky gets much
## less weight).
.asvf_level <- function(level) {
  list(poly_level = c(4, 8)[level], min_weight = c(0.4, 0.1)[level])
}

## Run the fused kernel over one tile (which already carries its overlap
## border). The reflect padding added here only affects the outermost band,
## which is either cropped away by the tiling engine or - at the true raster
## edge - is exactly the extension RVT uses.
##
## `negate`: negative openness is positive openness computed on the inverted
## terrain (RVT: `dem_arr = -1 * dem_arr` before the same horizon search) -
## depressions read as openness the way ridges do on the real terrain. No
## separate kernel needed, just flip the sign of what goes in.
.horizon_tile <- function(tile, xres, yres, num_directions, radius_max,
                           noise_removal, want_svf, want_opns, threads,
                           negate = FALSE, want_asvf = FALSE,
                           dir_weight = numeric(0)) {
  if (negate) tile <- -tile
  off <- .direction_offsets(num_directions, radius_max,
                             .radius_min(radius_max, noise_removal), xres, yres)
  pad <- as.integer(radius_max + 1L)
  padded <- .pad_reflect(tile, pad)
  out <- horizon_svf_opns(padded, pad, nrow(tile), ncol(tile),
                           off$dx, off$dy, off$dist, off$starts, off$ends,
                           want_svf, want_opns, want_asvf, dir_weight, threads)
  list(svf = out$svf, opns = out$opns, asvf = out$asvf)
}

## Multi-resolution variant of .horizon_tile(): level 0 is the tile, the rest
## are windows onto coarser copies of the DEM supplied by .process_tiled().
.horizon_pyramid_tile <- function(tile, ctx, xres, yres, offs, plan,
                                   want_svf, want_opns, want_asvf, dir_weight,
                                   threads, negate = FALSE) {
  if (negate) {
    tile <- -tile
    ctx$aux <- lapply(ctx$aux, function(a) { a$m <- -a$m; a })
  }
  pad <- as.integer(plan[[1]]$rmax + 1L)
  padded <- .pad_reflect(tile, pad)
  out <- horizon_pyramid(padded, pad, nrow(tile), ncol(tile),
                          as.integer(ctx$x0), as.integer(ctx$y0),
                          lapply(ctx$aux, `[[`, "m"),
                          vapply(ctx$aux, `[[`, integer(1), "fac"),
                          vapply(ctx$aux, `[[`, integer(1), "cx0"),
                          vapply(ctx$aux, `[[`, integer(1), "cy0"),
                          lapply(offs, `[[`, "dx"), lapply(offs, `[[`, "dy"),
                          lapply(offs, `[[`, "dist"), lapply(offs, `[[`, "starts"),
                          lapply(offs, `[[`, "ends"),
                          want_svf, want_opns, want_asvf, dir_weight, threads)
  list(svf = out$svf, opns = out$opns, asvf = out$asvf)
}

## `negate` for negative openness searches the terrain turned upside down, so
## every level is negated per tile. The coarse levels must then be built with
## `min` rather than `max`: the horizon wanted is max(-z), which is -min(z).
## Coarsening with max and negating would give -max(z) - the wrong extreme, and
## a level that hides obstructions instead of preserving them.
.horizon_run <- function(dem, out_path, want_svf, want_opns,
                          num_directions, reach, noise_removal,
                          tile_size, threads, overwrite, progress,
                          negate = FALSE, want_asvf = FALSE,
                          dir_weight = numeric(0),
                          pyramid_px = 100, pyramid_factor = 4) {
  if (!overwrite && file.exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  nm <- if (want_asvf) "asvf" else if (want_svf) "svf" else "opns"

  plan <- .pyramid_plan(info$xres, reach, pyramid_px, pyramid_factor)
  radius_max <- plan[[1]]$rmax

  # A reach that level 0 already covers needs no coarse levels at all - it is
  # exactly the ordinary search, so use it and keep one code path tested.
  if (length(plan) > 1L) {
    # noise_removal raises the inner edge of the full-resolution level; the
    # coarse levels start where it stops regardless.
    plan[[1]]$rmin <- .radius_min(radius_max, noise_removal)
    paths <- .pyramid_build(dem, plan, threads,
                             method = if (negate) "min" else "max")
    offs <- .pyramid_offsets(plan, num_directions, info$xres, info$yres)
    overlap <- as.integer(radius_max + 1L)
    if (is.null(tile_size))
      tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

    .process_tiled(dem, stats::setNames(list(out_path), nm), overlap, tile_size,
                    function(tile, xres, yres, ctx)
                      .horizon_pyramid_tile(tile, ctx, xres, yres, offs, plan,
                                             want_svf, want_opns, want_asvf,
                                             dir_weight, threads, negate),
                    progress = progress, threads = threads,
                    aux = .pyramid_aux(paths, plan))
    return(invisible(out_path))
  }

  overlap <- as.integer(radius_max + 1L)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

  .process_tiled(dem, stats::setNames(list(out_path), nm), overlap, tile_size,
                  function(tile, xres, yres)
                    .horizon_tile(tile, xres, yres, num_directions, radius_max,
                                   noise_removal, want_svf, want_opns, threads,
                                   negate, want_asvf, dir_weight),
                  progress = progress, threads = threads)
  invisible(out_path)
}

#' Sky-view factor
#'
#' How much of the sky each cell can see, from 0 (completely enclosed - the
#' bottom of a narrow pit or gorge) to 1 (completely open - flat ground or a
#' summit). Flat ground gives exactly 1.
#'
#' It looks like the terrain lit by a uniformly bright overcast sky. There is
#' no sun position, so nothing disappears into cast shadow and no feature is
#' favoured just because it happens to face the light - which makes it a good
#' general-purpose alternative to hillshade when features may run in any
#' direction. Hollows, ditches and the insides of banks come out dark; ridges,
#' mounds and open ground come out bright.
#'
#' @section How it works:
#' From each cell the terrain is scanned outwards along `num_directions`
#' evenly spaced compass directions (16 by default, i.e. every 22.5 degrees):
#'
#' * Along each direction, samples are taken out to `reach`,
#'   finely enough (three samples per pixel) that no pixel the ray crosses is
#'   skipped; duplicates are then dropped.
#' * For each sample, the angle from the cell up to that sample is computed
#'   from the elevation difference over the true ground distance. Pixel size
#'   is taken into account, so non-square pixels are handled correctly - but
#'   it also means **the elevation units must match the map units** (metres
#'   over metres) for the angles to be meaningful.
#' * The steepest angle found along a direction is that direction's horizon.
#'   Horizons below the horizontal are treated as horizontal, since a cell on
#'   a summit sees a full half-dome and no more.
#' * The result is the mean of `1 - sin(horizon angle)` over all directions.
#'
#' Everything runs at the raster's native resolution, tile by tile, with
#' tiles overlapping by the search radius so the result never depends on how
#' the raster was split up.
#'
#' @section Tuning:
#' * `reach` is the main control: it sets the size of the features the
#'   result responds to. It is in **map units** (metres, normally), so the
#'   same number means the same ground distance on any DEM. Short reaches
#'   emphasise fine local texture; long ones bring in shading from more
#'   distant landforms and take proportionally longer.
#' * `num_directions` trades smoothness for time (cost is roughly linear).
#'   16 is the usual choice; 32 reduces the faint directional banding that
#'   can appear in smooth terrain.
#' * `noise_removal` skips the innermost part of each ray. Raise it on
#'   noisy DEMs - vegetation or point-cloud speckle right next to a cell can
#'   otherwise set that cell's horizon on its own. Kokalj and Hesse note this
#'   markedly improves visibility of archaeological features where the data
#'   suffer from flight-strip misalignment or were gridded finer than the
#'   point density really supports.
#'
#' @section Recommended settings:
#' Kokalj and Hesse (2017) suggest, for archaeological work on any terrain, a
#' **10 m** search radius with **16 directions** - which is the default here,
#' `reach = 10`, whatever the resolution of the DEM. Much larger radii belong
#' to other disciplines - they cite 10 km for meteorological work - and only
#' generalise archaeological detail away.
#'
#' For display they suggest a linear stretch of 0.65-1.0 on varied terrain and
#' 0.9-1.0 on very flat ground, where the useful values are otherwise squeezed
#' into a narrow band just below 1.
#'
#' It works best on moderate to steep topography, where it brings out
#' depressions and features lying on slopes. On flat or very gentle ground it
#' is largely limited to negative features - pits, ditches, quarries, dolines -
#' and becomes very sensitive to DEM noise, which is where `noise_removal`
#' earns its keep.
#'
#' @references
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Prostor, kraj, čas 14. Ljubljana:
#' Založba ZRC. \doi{10.3986/9789612549848}
#'
#' Zakšek, K., Oštir, K. and Kokalj, Ž. (2011) Sky-View Factor as a Relief
#' Visualization Technique. *Remote Sensing* 3, 398-415.
#'
#' @param dem path to a DEM, or a vector of paths forming a mosaic (see
#'   [rvt_mosaic()]); mosaics are read across file boundaries, so tiles do not
#'   produce seams
#' @param out_path output GeoTIFF path (default: a temp file)
#' @param num_directions number of compass directions scanned (default 16)
#' @param reach how far to look, in **map units** (metres, normally;
#'   default 10). Scanned at full resolution for the first `pyramid_px`
#'   cells (100 by default, so 100 m on a 1 m DEM); past that the search
#'   reads from progressively coarser copies of the DEM, which is what keeps
#'   a long reach affordable - see [rvt_reach]
#' @param noise_removal 0-3; raises the minimum search radius to ignore
#'   near-cell noise (0, 10, 20 or 40 percent of `reach`)
#' @param reach look this far in **map units** (metres, normally) using a
#'   multi-resolution search, instead of `radius_max` pixels at full
#'   resolution. Cost then grows with the logarithm of the distance rather
#'   than in proportion to it, which is what makes a reach of hundreds of
#'   metres affordable on fine data. `radius_max` and `noise_removal` are
#'   ignored when this is set. See [rvt_reach].
#' @param pyramid_px,pyramid_factor shape of that multi-resolution search:
#'   pixels scanned per level (default 100, the accuracy knob) and how much
#'   coarser each level is than the last (default 4). Only used with `reach`.
#' @param tile_size tile edge in pixels; NULL picks a size targeting roughly
#'   256 MB per working matrix
#' @param threads C++ threads (see [rvt_threads()])
#' @param overwrite recompute even if `out_path` exists
#' @param progress report per-tile progress
#' @return `out_path`, invisibly. `dem` is the first argument, so this is
#'   pipe-friendly: `dem |> rvt_svf("svf.tif")`.
#' @seealso [rvt_openness()] for a related measure of the same horizon that
#'   spreads convex and concave terrain over a wider range.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_svf(dem)
#' @export
rvt_svf <- function(dem, out_path = tempfile(fileext = ".tif"),
                num_directions = 16, reach = 10, noise_removal = 0,
                pyramid_px = 100, pyramid_factor = 4,
                tile_size = NULL, threads = rvt_threads(),
                overwrite = FALSE, progress = FALSE) {
  .horizon_run(dem, out_path, TRUE, FALSE, num_directions, reach,
                noise_removal, tile_size, threads, overwrite, progress,
                pyramid_px = pyramid_px, pyramid_factor = pyramid_factor)
}

#' Anisotropic sky-view factor
#'
#' Sky-view factor computed under a sky that is brighter on one side, instead
#' of uniformly bright all round. Still 0-1, and still 1 on flat ground, but
#' terrain now casts a soft directional shading: slopes facing the bright
#' side come out lighter, slopes turned away from it darker.
#'
#' The point is to recover some of the readability of a hillshade - the eye is
#' good at reading relief lit from one side - without a hillshade's main flaw,
#' which is that features running parallel to the light almost disappear.
#' Nothing is ever fully in shadow here, so a feature aligned with
#' `main_direction` is dimmed rather than lost.
#'
#' @section How it works:
#' The horizon search is exactly the one in [rvt_svf()]. The only difference
#' is in the final averaging: instead of every direction counting equally,
#' each gets a weight running from 1 in `main_direction` down to `min_weight`
#' at the opposite azimuth, following `cos(half the angle from
#' main_direction)` raised to a power. The weighted contributions are divided
#' by the total weight, which keeps the result on the same 0-1 scale as plain
#' sky-view factor.
#'
#' `level` picks how sharply the weight falls away: level 1 is gentle
#' (opposite sky still counts 0.4), level 2 is strong (opposite sky counts
#' 0.1), which gives a more hillshade-like image.
#'
#' @section Tuning:
#' * `main_direction` is a compass azimuth in degrees - 0 north, 90 east, 180
#'   south, 270 west. The 315 default is north-west, the conventional
#'   illumination direction for terrain images, and matches
#'   [rvt_hillshade()]'s default so the two can be compared or blended.
#'   Rotate it to bring out features running in a particular direction: relief
#'   perpendicular to the light shows up most.
#' * `level` 1 keeps the result close to plain sky-view factor with a hint of
#'   direction; 2 pushes the directional effect much harder.
#' * `reach`, `num_directions` and `noise_removal` behave exactly as in
#'   [rvt_svf()]. Raising `num_directions` matters a little more here than for
#'   plain sky-view factor, since the weights vary between directions and a
#'   coarse set samples that variation coarsely.
#'
#' @section Recommended settings:
#' Radius and directions as for [rvt_svf()] - Kokalj and Hesse (2017) use
#' 10 m and 16 directions for archaeological work, and the same 0.65-1.0 /
#' 0.9-1.0 display stretches apply.
#'
#' Their assessment is that the anisotropy "brings back some of the
#' plasticity of hill shading and gives better details on very flat areas",
#' which is a fair summary of when to prefer this over plain sky-view factor:
#' reach for it when the isotropic image reads as flat and hard to interpret,
#' particularly on level ground. Both are described as very good general-
#' purpose choices, because they show small features whatever their
#' orientation and shape, on most terrain.
#'
#' @references
#' Zakšek, K., Oštir, K., Pehani, P., Kokalj, Ž. and Polert, E. (2012) Hill
#' Shading Based on Anisotropic Diffuse Illumination. In *Symposium GIS
#' Ostrava 2012*, 1-10. Ostrava: Technical University of Ostrava.
#'
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#'
#' @inheritParams rvt_svf
#' @param main_direction compass azimuth, in degrees, of the brightest part
#'   of the sky (default 315, i.e. north-west)
#' @param level strength of the anisotropy: 1 gentle (default), 2 strong
#' @return `out_path`, invisibly
#' @seealso [rvt_svf()] for the uniform-sky version, [rvt_hillshade()] for a
#'   conventional directional shading.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_asvf(dem)
#' @export
rvt_asvf <- function(dem, out_path = tempfile(fileext = ".tif"),
                num_directions = 16, reach = 10, noise_removal = 0,
                main_direction = 315, level = 1,
                pyramid_px = 100, pyramid_factor = 4,
                tile_size = NULL, threads = rvt_threads(),
                overwrite = FALSE, progress = FALSE) {
  if (!level %in% c(1, 2))
    stop("`level` must be 1 (gentle anisotropy) or 2 (strong)", call. = FALSE)
  a <- .asvf_level(level)
  w <- .asvf_weights(num_directions, main_direction, a$poly_level, a$min_weight)
  .horizon_run(dem, out_path, FALSE, FALSE, num_directions, reach,
                noise_removal, tile_size, threads, overwrite, progress,
                negate = FALSE, want_asvf = TRUE, dir_weight = w,
                pyramid_px = pyramid_px, pyramid_factor = pyramid_factor)
}

#' Positive openness
#'
#' The average angle from straight up (the zenith) down to the horizon, in
#' degrees. Flat ground gives exactly 90 degrees. Above 90 means the cell
#' stands proud of its surroundings - ridges, mounds, bank tops, wall lines;
#' below 90 means it sits enclosed - ditches, hollows, the inside of a bank.
#'
#' This measures the same horizon as [rvt_svf()] but keeps more of its range.
#' Sky-view factor compresses everything above the horizontal into a narrow
#' band just below 1, whereas openness spreads convex, flat and concave
#' terrain out evenly around 90 degrees - so the edges of low earthworks tend
#' to come out crisper. Values are not restricted to 0-90: around 107 degrees
#' is reached on strongly convex terrain in the bundled sample tile.
#'
#' @section How it works:
#' The horizon search is identical to [rvt_svf()] - same directions, same
#' sampling along each ray, same use of true ground distance - with two
#' differences in how the horizon angles are combined:
#'
#' * Horizons below the horizontal are **not** clamped. A cell looking down
#'   into falling ground keeps its negative angle, which is what lets convex
#'   terrain rise above 90 degrees.
#' * The result is `90 - mean(horizon angle)` in degrees, over all directions.
#'
#' Like sky-view factor it runs at native resolution, tile by tile, with
#' tiles overlapping by the search radius. Elevation units must match the map
#' units for the angles to mean anything.
#'
#' @section Tuning:
#' The parameters behave exactly as in [rvt_svf()] - `reach` sets the size of
#' features the result responds to, in map units, `num_directions` trades
#' smoothness against time, and `noise_removal` protects against speckle
#' immediately next to a cell.
#'
#' A common pairing is to compute this alongside [rvt_openness_negative()]:
#' the positive image favours convex features, the negative one concave, and
#' comparing or differencing the two separates banks from ditches more
#' clearly than either does alone.
#'
#' @section Recommended settings:
#' As for [rvt_svf()]: Kokalj and Hesse (2017) use a **10 m** radius
#' (`reach = 10`) with 16 directions, and suggest
#' displaying positive openness with a linear stretch of 65-95 degrees on
#' varied terrain, 85-91 on very flat ground, 55-95 on steep or complex
#' ground.
#'
#' Because openness ignores the mathematical horizon and takes in the whole
#' sphere, the image comes out notably *flat* - the broad shape of the
#' landscape largely disappears, rather as if the trend had been removed.
#' That costs you the intuitive read of the terrain that sky-view factor
#' keeps, and makes interpretation a little harder, but it buys two things:
#' the result is not washed out by gentle or steep slopes and so works across
#' varied topography, and a given feature looks much the same wherever it sits
#' - on the flat or partway up a hillside. Doneus (2013) makes the point that
#' this consistency of "signature" is what makes openness well suited to
#' automated feature detection, not just to looking.
#'
#' @references
#' Yokoyama, R., Shirasawa, M. and Pike, R. J. (2002) Visualizing Topography
#' by Openness: A New Application of Image Processing to Digital Elevation
#' Models. *Photogrammetric Engineering and Remote Sensing* 68, 251-266.
#'
#' Doneus, M. (2013) Openness as Visualization Technique for Interpretative
#' Mapping of Airborne LiDAR Derived Digital Terrain Models. *Remote Sensing*
#' 5, 6427-6442.
#'
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#'
#' @inheritParams rvt_svf
#' @return `out_path`, invisibly
#' @seealso [rvt_openness_negative()] for the concave counterpart,
#'   [rvt_svf()] for the 0-1 sky-fraction version of the same horizon.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_openness(dem)
#' @export
rvt_openness <- function(dem, out_path = tempfile(fileext = ".tif"),
                     num_directions = 16, reach = 10, noise_removal = 0,
                     pyramid_px = 100, pyramid_factor = 4,
                     tile_size = NULL, threads = rvt_threads(),
                     overwrite = FALSE, progress = FALSE) {
  .horizon_run(dem, out_path, FALSE, TRUE, num_directions, reach,
                noise_removal, tile_size, threads, overwrite, progress,
                pyramid_px = pyramid_px, pyramid_factor = pyramid_factor)
}

#' Negative openness
#'
#' The same measurement as [rvt_openness()], made on the terrain turned upside
#' down. Concave features become the prominent ones: ditches, pits, hollow
#' ways and quarry scoops read high, where positive openness favours banks,
#' mounds and ridges.
#'
#' Note that only the *shape* preference is inverted, not the numbers. Flat
#' ground still gives exactly 90 degrees, and a cell in a pit reads *above*
#' 90 here - because on the inverted surface that pit is a mound. Values are
#' in degrees and, as with positive openness, are not restricted to 0-90.
#'
#' @section How it works:
#' Elevations are negated, and then exactly the horizon search described in
#' [rvt_openness()] runs on the result - same directions, same sampling, same
#' `90 - mean(horizon angle)`. Nothing else differs, so the two are directly
#' comparable and can be differenced.
#'
#' @section Tuning:
#' Identical to [rvt_openness()]. Use the same `reach` for both if you
#' intend to compare or difference them - a mismatch there makes the two
#' images respond to different feature sizes and the comparison stops meaning
#' much.
#'
#' @section Recommended settings:
#' Worth being clear about, since the name invites the assumption: **negative
#' openness is not the inverse of positive openness.** It carries genuinely
#' different information rather than the same image with the sign flipped.
#' Where positive openness picks out the convexities - the ridge *between* two
#' hollow ways, the raised rim of a crater - negative openness picks out the
#' deepest parts of the concavities: the hollow way itself, the floor of a
#' gorge, the foot of a cliff. Computing both and reading them together is the
#' point.
#'
#' Kokalj and Hesse (2017) recommend displaying it with an **inverted**
#' greyscale, dark for high values, so that concave features read dark in both
#' images and the pair stay visually consistent. Their suggested stretch is
#' 60-95 degrees on varied terrain, 75-95 on very flat ground and 45-95 on
#' steep or complex ground. Radius and directions as for [rvt_svf()]: 10 m,
#' 16 directions.
#'
#' @references
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#'
#' @inheritParams rvt_svf
#' @return `out_path`, invisibly
#' @seealso [rvt_openness()] for the convex counterpart.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_openness_negative(dem)
#' @export
rvt_openness_negative <- function(dem, out_path = tempfile(fileext = ".tif"),
                     num_directions = 16, reach = 10, noise_removal = 0,
                     pyramid_px = 100, pyramid_factor = 4,
                     tile_size = NULL, threads = rvt_threads(),
                     overwrite = FALSE, progress = FALSE) {
  .horizon_run(dem, out_path, FALSE, TRUE, num_directions, reach,
                noise_removal, tile_size, threads, overwrite, progress,
                negate = TRUE,
                pyramid_px = pyramid_px, pyramid_factor = pyramid_factor)
}
