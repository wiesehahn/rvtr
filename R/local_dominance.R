## Local dominance (Hesse's LiVT algorithm, as used by RVT).

## Grid of (angle, distance) samples, flattened to integer pixel offsets - no
## deduplication and no direction-by-direction structure like the horizon
## kernel's, because local dominance sums *every* (angle, distance) sample
## independently rather than taking a max per direction. Multiple samples can
## round to the same integer offset (small distances, fine angular_res); RVT
## doesn't deduplicate those either, so neither do we - it would change the
## weighting.
##
## Angle 0 is arbitrary and the offset formula doesn't negate the y term the
## way .direction_offsets() does for the horizon kernel: the angle sweep is a
## full 0-360 circle at a fixed step, which is symmetric under negating either
## component, so the aggregate (summed) result is identical either way.
.local_dominance_offsets <- function(min_rad, max_rad, rad_inc, angular_res) {
  distances <- seq(min_rad, max_rad, by = rad_inc)
  angles <- seq(0, 360 - angular_res, by = angular_res) * pi / 180
  grid <- expand.grid(angle = angles, distance = distances)
  list(dx = as.integer(round(cos(grid$angle) * grid$distance)),
       dy = as.integer(round(sin(grid$angle) * grid$distance)),
       dist = grid$distance,
       weight = 2 * grid$distance + rad_inc,
       distances = distances, n_ang = length(angles))
}

.local_dominance_tile <- function(tile, min_rad, max_rad, rad_inc, angular_res,
                                   observer_height, threads) {
  off <- .local_dominance_offsets(min_rad, max_rad, rad_inc, angular_res)
  pad <- as.integer(max_rad)
  padded <- .pad_edge(tile, pad)
  out <- local_dominance_kernel(padded, pad, nrow(tile), ncol(tile),
                                 off$dx, off$dy, off$dist, off$weight,
                                 observer_height, threads)
  norma <- sum((observer_height / off$distances) * (2 * off$distances + rad_inc)) * off$n_ang
  list(ld = out / norma)
}

#' Local dominance
#'
#' How much the terrain around each cell falls away below an observer standing
#' on it. High where you would look down on your surroundings - mounds, ridge
#' crests, plateau edges, the tops of banks - and low in hollows, ditches and
#' at the foot of slopes.
#'
#' The result is scaled so that **flat ground gives exactly 1**: above 1 means
#' more dominant than flat, below 1 means overlooked by the surroundings. On
#' the bundled sample tile values run from roughly 0.2 to 2.6.
#'
#' Where [rvt_openness()] only cares about the highest thing on each horizon,
#' this samples the whole surrounding ring and weights every point by how far
#' below eye level it sits. Broad, gentle rises therefore register here even
#' though nothing about them breaks a horizon - which is what makes it good at
#' picking out low, spread-out features such as plateau edges and slight
#' mounds.
#'
#' @section How it works:
#' Around each cell, terrain is sampled on a ring rather than along rays to
#' the horizon:
#'
#' * Sample points sit between `min_rad` and `max_rad` out, stepping
#'   by `rad_inc`, and every `angular_res` degrees around the full circle -
#'   by default 11 distances times 24 directions, so 264 samples per cell.
#'   Each is rounded to the nearest pixel.
#' * The observer's eye is placed at the cell's elevation plus
#'   `observer_height`. Samples at or above eye level contribute nothing; each
#'   one below adds its drop below the eye divided by its distance.
#' * Each contribution is weighted in proportion to the area of the ring it
#'   came from, so distant rings - which contain more ground - are not
#'   drowned out by the handful of pixels close in.
#' * The total is divided by what perfectly flat ground would produce, which
#'   is what puts flat terrain at exactly 1.
#'
#' Distances here are in map units, as everywhere else, so the same settings
#' mean the same ground distances on any DEM. Computation is at native
#' resolution, tile by tile, with tiles overlapping by `max_rad`.
#'
#' @section Tuning:
#' * `min_rad` and `max_rad` set the band of feature sizes the result responds
#'   to, and they are in **map units** (metres, normally): the defaults cover
#'   10-20 m on any DEM. Widen the band for larger landforms.
#' * Starting at `min_rad` rather than at the cell itself is deliberate - it
#'   keeps immediate micro-relief and DEM noise from dominating. Lower it to
#'   pick up smaller features, raise it to ignore more local roughness.
#' * `angular_res` controls how many directions are sampled (360 divided by
#'   it). Smaller values give a smoother, less spoke-like result and cost
#'   proportionally more time.
#' * `rad_inc` controls sampling density along each direction. Larger steps
#'   are faster and coarser.
#' * `observer_height` sets how sensitive the result is: a low observer
#'   exaggerates small height differences, a tall one flattens the contrast
#'   because everything already lies well below the eye. The 1.7 m default is
#'   roughly human eye level.
#'
#' @section Recommended settings:
#' Kokalj and Hesse (2017) recommend a **10-20 m** search band with a 1.7 m
#' observer on flat to moderate terrain, narrowing to about 10 m on steep or
#' complex ground. Those are metres, which is what these parameters take, so
#' the defaults already match their advice on any DEM. For display they
#' suggest a linear stretch of roughly 0.5-1.8,
#' opening out to 0.5-3.0 on very flat ground.
#'
#' It is singled out as one of the better choices for **very subtle positive
#' relief** - former field boundaries, ploughed-out burial mounds - while also
#' doing well on depressions such as dolines, mining traces and hollow ways.
#'
#' One property to keep in mind when picking a stretch: unlike openness, this
#' keeps a limited impression of the overall landscape, because a cell part
#' way down a slope genuinely does look down on more ground than one on the
#' flat. Steep country therefore sits at systematically higher values and
#' generally wants a different stretch from gentle country - which is also
#' why flat ground reading exactly 1 is a useful anchor.
#'
#' @references
#' Hesse, R. (2016) Visualisierung hochauflösender digitaler
#' Geländemodelle mit LiVT. In *Computeranwendungen und quantitative Methoden
#' in der Archäologie. 4. Workshop der AG CAA 2013*, edited by U. Lieberwirth
#' and I. Herzog, 109-128. Berlin: Topoi.
#'
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#'
#' @inheritParams rvt_svf
#' @param min_rad,max_rad inner/outer radius of the sampled ring, in **map
#'   units** (metres, normally)
#'   (defaults 10 and 20)
#' @param rad_inc distance between samples along each ray, in map units;
#'   anything finer than one cell means every cell (default 1)
#' @param angular_res angle between sampled directions, in degrees
#'   (default 15, giving 24 directions)
#' @param observer_height height of the observer above the terrain, in the
#'   DEM's elevation units, typically metres (default 1.7)
#' @return `out_path`, invisibly
#' @seealso [rvt_openness()], which responds to similar terrain but only
#'   through the horizon.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_local_dominance(dem)
#' @export
rvt_local_dominance <- function(dem, out_path = fs::file_temp(ext = "tif"),
                                 min_rad = 10, max_rad = 20, rad_inc = 1,
                                 angular_res = 15, observer_height = 1.7,
                                 tile_size = NULL, threads = rvt_threads(),
                                 overwrite = FALSE, progress = FALSE) {
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  max_px <- .cells(max_rad, info$xres, "max_rad", "outer")
  min_px <- .cells(min_rad, info$xres, "min_rad", "inner")
  inc_px <- .cells(rad_inc, info$xres, "rad_inc", "step")
  if (min_px > max_px)
    stop("`min_rad` must not exceed `max_rad`", call. = FALSE)

  overlap <- as.integer(max_px)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

  .process_tiled(dem, list(ld = out_path), overlap, tile_size,
                  function(tile, xres, yres)
                    .local_dominance_tile(tile, min_px, max_px, inc_px,
                                           angular_res, observer_height, threads),
                  progress = progress, threads = threads)
  invisible(out_path)
}
