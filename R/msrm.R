## Multi-scale relief model.

## Sequence of box-filter radii (in pixels): `ndx^scaling_factor` for ndx
## running from the smallest feature size up to the largest, both converted
## from map units to pixels via `resolution`. feature_min is clamped up to
## resolution - a feature smaller than one pixel isn't representable.
.msrm_radii <- function(feature_min, feature_max, scaling_factor, resolution) {
  feature_min <- max(feature_min, resolution)
  i <- floor(((feature_min - resolution) / (2 * resolution)) ^ (1 / scaling_factor))
  n <- ceiling(((feature_max - resolution) / (2 * resolution)) ^ (1 / scaling_factor))
  as.integer(round((i:n) ^ scaling_factor))
}

.msrm_tile <- function(tile, radii, threads) {
  pad <- as.integer(max(radii))
  padded <- .pad_edge(tile, pad)
  out <- msrm_kernel(padded, pad, nrow(tile), ncol(tile), radii, threads)
  list(msrm = out)
}

## Simple local relief model is the same kernel with just two radii: a box
## mean of radius 0 is the cell itself, so the single "difference between
## consecutive smoothed surfaces" it computes is exactly
## `dem - mean_filter(dem, radius)`.
.slrm_tile <- function(tile, radius, threads) {
  pad <- as.integer(radius)
  padded <- .pad_edge(tile, pad)
  out <- msrm_kernel(padded, pad, nrow(tile), ncol(tile),
                      c(0L, as.integer(radius)), threads)
  list(slrm = out)
}

#' Multi-scale relief model (MSRM)
#'
#' Strips away both the broad shape of the landscape and fine-grained noise,
#' leaving only relief within a chosen size range.
#'
#' Values are in the DEM's elevation units and read directly as height above
#' or below the local surface: positive for anything raised (banks, mounds,
#' walls), negative for anything cut in (ditches, hollow ways, pits). A value
#' of 0.3 on a metric DEM means the cell sits about 0.3 m proud of its
#' surroundings at the scales you selected - which makes this unusually easy
#' to interpret compared with the angle- and ratio-based metrics here.
#'
#' Because the regional trend is removed, a low bank on a steep hillside shows
#' up just as clearly as one on level ground - the situation where hillshade
#' and sky-view factor struggle most. **Flat ground gives 0, and so does a
#' uniform slope**: only departures from the local trend survive.
#'
#' @section How it works:
#' The DEM is smoothed repeatedly, by averaging over square windows of
#' increasing size, and it is the *differences between* those smoothed
#' versions that form the result:
#'
#' * The window radii come from `feature_min` and `feature_max`, converted
#'   from map units to pixels using the raster's own resolution, and spaced
#'   apart as `n^scaling_factor`. At the defaults on a 1 m DEM that gives
#'   radii of 0, 1, 4, 9 and 16 pixels.
#' * Subtracting each smoothed surface from the previous, less-smoothed one
#'   leaves exactly the detail that step of smoothing removed - the relief at
#'   that particular scale.
#' * Those differences are averaged into a single band-pass image, keeping
#'   what lies between the smallest and largest scale and discarding both
#'   finer noise and broader landscape form.
#'
#' The window averages are computed from a summed-area table, so their cost
#' does not grow with window size - a large `feature_max` is very nearly as
#' cheap as a small one. NoData cells are excluded from the averages rather
#' than counted as zero. Computation is at native resolution, tile by tile,
#' with tiles overlapping by the largest radius.
#'
#' @section Tuning:
#' * `feature_min` is the smallest thing kept. Leaving it at 0 (it is clamped
#'   up to one pixel) keeps the sharpest detail; raise it to a few times the
#'   pixel size on noisy lidar or vegetation-affected DEMs to suppress
#'   speckle.
#' * `feature_max` sets how much landscape form is removed - features much
#'   larger than it fade out. Lower it to flatten strong topography harder
#'   and concentrate on small features; raise it to keep broader ones.
#' * `scaling_factor` controls how many scales are sampled between the two.
#'   1 samples every radius: the most scales, smoothest response, best
#'   sensitivity to intermediate feature sizes, slightly slower. The default
#'   2 uses fewer, more widely spaced scales, which gives more contrast.
#'   Higher values are punchier still but can step over intermediate sizes.
#' * Note `feature_min`/`feature_max` are in **map units, not pixels** (the
#'   opposite of [rvt_local_dominance()]), so the same settings mean the same
#'   ground distances regardless of resolution.
#'
#' @inheritParams rvt_svf
#' @param feature_min,feature_max smallest/largest feature size to keep, in
#'   the DEM's horizontal map units, typically metres; `feature_min` is
#'   clamped up to the pixel resolution if smaller (defaults 0 and 20)
#' @param scaling_factor positive integer spacing the sampled radii
#'   (`radius = n^scaling_factor`); 1 samples every radius, higher values
#'   sample fewer, more widely spaced ones (default 2)
#' @return `out_path`, invisibly
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_msrm(dem)
#' @export
rvt_msrm <- function(dem, out_path = tempfile(fileext = ".tif"),
                      feature_min = 0, feature_max = 20, scaling_factor = 2,
                      tile_size = NULL, threads = rvt_threads(),
                      overwrite = FALSE, progress = FALSE) {
  if (!overwrite && file.exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  radii <- .msrm_radii(feature_min, feature_max, scaling_factor, info$xres)
  if (length(unique(radii)) < 2)
    stop("feature_min/feature_max/scaling_factor produce fewer than two ",
         "distinct radii at this resolution; widen the range", call. = FALSE)

  overlap <- as.integer(max(radii))
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

  .process_tiled(dem, list(msrm = out_path), overlap, tile_size,
                  function(tile, xres, yres) .msrm_tile(tile, radii, threads),
                  progress = progress, threads = threads)
  invisible(out_path)
}

#' Simple local relief model (SLRM)
#'
#' Height of each cell above or below the smoothed terrain around it, in the
#' DEM's elevation units. Positive where the ground is locally raised (banks,
#' mounds, walls), negative where it is cut in (ditches, hollow ways, pits),
#' and 0 where it follows the local trend.
#'
#' This is the simplest of the trend-removal visualizations and still one of
#' the most used: subtracting a smoothed copy of the terrain from itself
#' strips out hillslopes and valleys, so low features read the same on a
#' slope as on level ground. [rvt_msrm()] does the same thing across a range
#' of scales at once and is more forgiving when you don't know the size of
#' what you are looking for; this one gives you a single, sharply defined
#' scale that you set directly.
#'
#' @section How it works:
#' The DEM is averaged over a square window of `radius` pixels, and that
#' smoothed surface is subtracted from the original. What remains is
#' everything smaller than roughly the window - anything broader has been
#' absorbed into the average and cancels out.
#'
#' The window average comes from a summed-area table, so cost does not grow
#' with the radius, and NoData cells are left out of the average rather than
#' counted as zero. Computation is at native resolution, tile by tile, with
#' tiles overlapping by `radius`.
#'
#' @section Tuning:
#' * `radius` is the one control, and it is in **pixels**: the default 20 is
#'   20 m on a 1 m DEM but only 5 m on a 0.25 m DEM, so scale it with your
#'   resolution. It should be comfortably larger than the features you want
#'   to keep - features approaching the window size get partly absorbed into
#'   the trend and fade out.
#' * Too small a radius flattens everything and leaves only noise; too large
#'   leaves the hillslopes in and defeats the point. Somewhere around two to
#'   three times the width of your target features is a reasonable starting
#'   point.
#' * A square averaging window leaves faint blocky artefacts around very
#'   sharp features. If that becomes distracting, [rvt_msrm()] averages
#'   several scales and largely avoids it.
#'
#' @section Recommended settings:
#' Kokalj and Hesse (2017) suggest a filter radius of about **10 m** generally,
#' 5 m on very flat ground and 25 m in steep or complex terrain - metres
#' again, against a `radius` in pixels, so multiply by four on a 0.25 m DEM.
#' For display they use a linear stretch of -1 to +1 m, tightening to -0.5 to
#' +0.5 m on very flat ground and opening to -2 to +2 m on steep ground; the
#' right choice really depends on how tall the features you are after are.
#'
#' Trend removal is picked out as one of the most useful approaches for very
#' low relief on flat to moderate ground - former field boundaries, levelled
#' burial mounds. Compared with the closely related local relief model, this
#' simple form is faster and simpler, at the cost of somewhat less faithful
#' relative heights for the anomalies it reveals.
#'
#' @references
#' Hesse, R. (2010) LiDAR-derived Local Relief Models - a New Tool for
#' Archaeological Prospection. *Archaeological Prospection* 17, 67-72.
#'
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#'
#' @inheritParams rvt_svf
#' @param radius radius of the smoothing window, in *pixels* (default 20)
#' @return `out_path`, invisibly
#' @seealso [rvt_msrm()] for the multi-scale version.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_slrm(dem)
#' @export
rvt_slrm <- function(dem, out_path = tempfile(fileext = ".tif"),
                      radius = 20, tile_size = NULL, threads = rvt_threads(),
                      overwrite = FALSE, progress = FALSE) {
  radius <- as.integer(radius)
  if (is.na(radius) || radius < 1)
    stop("`radius` must be a positive number of pixels", call. = FALSE)
  if (!overwrite && file.exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, radius)

  .process_tiled(dem, list(slrm = out_path), radius, tile_size,
                  function(tile, xres, yres) .slrm_tile(tile, radius, threads),
                  progress = progress, threads = threads)
  invisible(out_path)
}
