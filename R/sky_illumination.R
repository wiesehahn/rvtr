## Sky illumination and cast shadow - both horizon searches, the first over
## every direction, the second along one.

## Compass azimuth of each of .direction_offsets()' directions. Those run
## counterclockwise from east, so direction d faces 90 - angle_d.
.direction_azimuths <- function(num_directions) {
  (90 - (360 / num_directions) * (0:(num_directions - 1))) %% 360
}

## Ray offsets for arbitrary compass azimuths, in the same flattened layout
## .direction_offsets() produces. Used where the direction is given rather
## than evenly spaced - a sun position, for instance.
.ray_offsets <- function(azimuth_deg, radius_max, radius_min, xres, yres,
                          oversample = 3) {
  radii <- seq(radius_min, radius_max, by = 1 / oversample)
  per_dir <- lapply(azimuth_deg, function(a) {
    ang <- (90 - a) * pi / 180
    dx <- round(cos(ang) * radii)
    dy <- round(-sin(ang) * radii)
    keep <- !duplicated(paste(dx, dy)) & !(dx == 0 & dy == 0)
    dx <- dx[keep]; dy <- dy[keep]
    ord <- order(sqrt(dx^2 + dy^2))
    list(dx = dx[ord], dy = dy[ord])
  })
  dx <- unlist(lapply(per_dir, `[[`, "dx"))
  dy <- unlist(lapply(per_dir, `[[`, "dy"))
  ends <- cumsum(vapply(per_dir, function(o) length(o$dx), integer(1)))
  list(dx = as.integer(dx), dy = as.integer(dy),
       dist = sqrt((dx * xres)^2 + (dy * yres)^2),
       starts = as.integer(c(0L, utils::head(ends, -1))),
       ends = as.integer(ends))
}

.sky_illumination_tile <- function(tile, ctx, xres, yres, offs, pad,
                                    num_directions, overcast, threads) {
  out <- do.call(sky_illumination_kernel,
                  c(list(.pad_reflect(tile, pad), pad,
                          nrow(tile), ncol(tile), xres, yres),
                    .level_args(ctx, offs),
                    list(.direction_azimuths(num_directions), overcast, threads)))
  list(sim = out)
}

#' Sky illumination
#'
#' How much diffuse daylight each cell receives from the sky as a whole -
#' accounting both for the way the ground is tilted and for the surrounding
#' terrain blocking part of the view. Scaled so flat, open ground reads 1;
#' below that means shaded by slope or surroundings.
#'
#' Unlike [rvt_hillshade()] there is no sun and no cast shadow, and unlike
#' [rvt_svf()] the sky is not counted uniformly: what a surface actually
#' receives depends on how squarely it faces each part of the sky, so a slope
#' tilted away from open sky is darker even with an identical horizon. The
#' result looks like an overcast day, which is a fair approximation of how
#' terrain reads to the eye in soft light.
#'
#' This is the most expensive metric here by a wide margin - see Tuning.
#'
#' @section How it works:
#' The horizon is found exactly as in [rvt_svf()] - a maximum slope along each
#' of `num_directions` rays out to `reach`. Each direction then
#' contributes light according to how much sky it leaves open above the
#' horizon *and* how that sky sits relative to the surface: the cell's slope
#' and aspect weight every direction, so sky behind the surface contributes
#' nothing. The `sky_model` sets how brightness varies with height in the sky
#' - `"uniform"` treats it as equally bright everywhere, `"overcast"` makes it
#' brighter overhead, which is the standard cloudy-sky distribution and the
#' default.
#'
#' Results are scaled by the value flat ground produces, so flat reads 1
#' regardless of `sky_model` or `num_directions`.
#'
#' @inheritSection rvt_daylight Direct sun versus diffuse sky
#'
#' @section Tuning:
#' * `reach` dominates the cost: the default 100 m means around 3500 terrain
#'   samples per cell, some 20 times what [rvt_svf()] does at its own default.
#'   Expect around a minute for a 4000 x 4000 raster. Shortening it to 30 m
#'   still captures the terrain that matters most for local shading; going
#'   *longer* is much cheaper than it used to be, since past one
#'   full-resolution scan - `pyramid_px` cells, 100 by default - the search
#'   coarsens as it goes (see [rvt_reach]).
#' * `num_directions` costs time in direct proportion. 32 is the default;
#'   16 halves the work at some loss of smoothness.
#' * `sky_model = "uniform"` is slightly cheaper and gives a flatter,
#'   more even image; `"overcast"` gives more contrast between open and
#'   enclosed ground.
#'
#' @section Recommended settings:
#' Kokalj and Hesse (2017) mention diffuse insolation only in passing, and
#' their verdict is worth repeating before you spend the runtime: it does
#' reveal archaeological features, but it needs considerably more computation
#' than the alternatives and the results come out **more generalised** than,
#' say, sky-view factor. If you want a diffuse-light view of the terrain and
#' are not specifically after the physical illumination quantity, [rvt_svf()]
#' gets you most of the way for a twentieth of the work, and [rvt_asvf()]
#' adds back a sense of direction.
#'
#' @references
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#'
#' @inheritParams rvt_svf
#' @param sky_model `"overcast"` (default) or `"uniform"` - how sky
#'   brightness varies with height above the horizon
#' @param num_directions number of directions searched (default 32)
#' @param reach how far to look, in **map units** (metres, normally;
#'   default 100). Scanned at full resolution for the first `pyramid_px`
#'   cells (100 by default, so 100 m on a 1 m DEM); past that the search
#'   reads from progressively coarser copies of the DEM, so a long reach
#'   stays affordable - see [rvt_reach]
#' @return `out_path`, invisibly
#' @seealso [rvt_daylight()] for direct sun over time rather than diffuse light
#'   at no particular time, and [rvt_svf()] for a much cheaper measure of the
#'   same horizon that ignores which way the ground faces.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' # small radius to keep the example quick
#' rvt_sky_illumination(dem, reach = 10, num_directions = 8)
#' @export
rvt_sky_illumination <- function(dem, out_path = fs::file_temp(ext = "tif"),
                                  sky_model = c("overcast", "uniform"),
                                  num_directions = 32, reach = 100,
                                  pyramid_px = 100, pyramid_factor = 4,
                                  pyramid_method = "max",
                                  tile_size = NULL, threads = rvt_threads(),
                                  overwrite = FALSE, progress = FALSE) {
  sky_model <- match.arg(sky_model)
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  py <- .pyramid_setup(dem, info, reach, pyramid_px, pyramid_factor, threads,
                        method = pyramid_method,
                        function(rmax, rmin, res)
                          .direction_offsets(num_directions, rmax, rmin, res, res))
  if (is.null(tile_size))
    tile_size <- .auto_tile_size(info$nx, info$ny, py$overlap)
  pad <- py$overlap

  .process_tiled(dem, list(sim = out_path), py$overlap, tile_size,
                  function(tile, xres, yres, ctx = .no_aux)
                    .sky_illumination_tile(tile, ctx, xres, yres, py$offs, pad,
                                            num_directions,
                                            sky_model == "overcast", threads),
                  progress = progress, threads = threads, aux = py$aux)
  invisible(out_path)
}

## A cast shadow is a one-direction horizon search, so it reuses the horizon
## kernel wholesale - including its multi-resolution levels, which matters here
## because a shadow-casting ridge is often a long way off.
.shadow_tile <- function(tile, ctx, xres, yres, offs, pad, sun_elevation,
                          threads) {
  h <- do.call(horizon_pyramid,
                c(list(.pad_reflect(tile, pad), pad, nrow(tile), ncol(tile)),
                  .level_args(ctx, offs),
                  list(FALSE, TRUE, FALSE, numeric(0), threads)))
  # the openness output of a one-direction search is 90 - horizon angle
  horizon <- 90 - h$opns
  list(shadow = (horizon < sun_elevation) * 1)
}

#' Cast shadow
#'
#' Which ground is lit and which lies in shadow for a given sun position: 1
#' where the sun reaches the cell, 0 where terrain between the cell and the
#' sun blocks it.
#'
#' This is genuine cast shadow, not the shading in [rvt_hillshade()]. Hillshade
#' only asks whether a surface is tilted towards the light, so it never puts
#' anything in the shadow of a hill in front of it; this traces the terrain
#' towards the sun and finds out. Low sun angles throw long shadows that pick
#' out very slight relief, which is why raking light is a standard trick for
#' spotting earthworks.
#'
#' @section How it works:
#' A single ray is traced from each cell towards `sun_azimuth`, out to
#' `reach`, keeping the steepest angle up to the terrain along it -
#' the horizon in the sun's direction. If that horizon stands higher than
#' `sun_elevation`, the sun is blocked and the cell is in shadow.
#'
#' Only one direction is searched, so this is far cheaper than the all-round
#' metrics: roughly one thirty-second of [rvt_sky_illumination()] at the same
#' radius. The output is 0 or 1, with NoData preserved.
#'
#' @inheritSection rvt_daylight Direct sun versus diffuse sky
#'
#' @section Tuning:
#' * `reach` sets how far a shadow can reach. A shadow cast by terrain
#'   further away than this is missed, so it needs to cover the tallest
#'   feature divided by the tangent of the sun elevation - a 20 m rise at 10
#'   degrees casts a shadow about 113 m long. Low sun therefore needs a large
#'   radius to be honest.
#' * `sun_elevation` is the main control on the picture: low values (5-15)
#'   throw long shadows that reveal very subtle relief, high values shadow
#'   only steep ground.
#' * `sun_azimuth` is a compass azimuth, as in [rvt_hillshade()].
#'
#' @inheritParams rvt_svf
#' @param sun_azimuth compass azimuth of the sun, in degrees (default 315)
#' @param sun_elevation height of the sun above the horizon, in degrees
#'   (default 35)
#' @param reach how far along the ray to look for blocking terrain, in
#'   **map units** (metres, normally; default 100). Scanned at full
#'   resolution for the first `pyramid_px` cells (100 by default, so 100 m on
#'   a 1 m DEM); past that the search reads from progressively coarser copies
#'   of the DEM, so a shadow cast from kilometres away stays affordable - see
#'   [rvt_reach]
#' @return `out_path`, invisibly - 1 lit, 0 shadowed
#' @seealso [rvt_daylight()], which integrates this over many sun positions
#'   instead of freezing one, and [rvt_hillshade()], which shades by surface
#'   orientation without casting shadows at all.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_shadow(dem, sun_elevation = 15, reach = 30)
#' @export
rvt_shadow <- function(dem, out_path = fs::file_temp(ext = "tif"),
                        sun_azimuth = 315, sun_elevation = 35, reach = 100,
                        pyramid_px = 100, pyramid_factor = 4,
                        pyramid_method = "max",
                        tile_size = NULL, threads = rvt_threads(),
                        overwrite = FALSE, progress = FALSE) {
  if (length(sun_azimuth) != 1L)
    stop("`sun_azimuth` must be a single direction", call. = FALSE)
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  py <- .pyramid_setup(dem, info, reach, pyramid_px, pyramid_factor, threads,
                        method = pyramid_method,
                        function(rmax, rmin, res)
                          .ray_offsets(sun_azimuth, rmax, rmin, res, res))
  if (is.null(tile_size))
    tile_size <- .auto_tile_size(info$nx, info$ny, py$overlap)
  pad <- py$overlap

  .process_tiled(dem, list(shadow = out_path), py$overlap, tile_size,
                  function(tile, xres, yres, ctx = .no_aux)
                    .shadow_tile(tile, ctx, xres, yres, py$offs, pad,
                                  sun_elevation, threads),
                  progress = progress, threads = threads, aux = py$aux)
  invisible(out_path)
}
