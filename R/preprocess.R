## Preparing a surface before measuring it: filling holes, removing noise.

#' Fill gaps in a raster
#'
#' Interpolates across NoData holes so the surface is continuous. Worth doing
#' before anything else, because a hole is not simply "no answer here": it
#' spreads. The derivative metrics ([rvt_slope()], [rvt_curvature()]) lose a
#' ring of cells around every hole, and the horizon metrics treat a hole as
#' *transparent* - light passes straight through ground that was never
#' measured - so a DSM with gaps in its canopy will read as brighter than the
#' real thing.
#'
#' @section How it works:
#' GDAL's `GDALFillNodata()` (via [gdalraster::fillNodata()]): for each cell
#' in a hole, the value is interpolated from the nearest valid cells in each
#' of four directions, weighted by inverse distance, followed by optional
#' smoothing passes over the interpolated values only.
#'
#' Cells that already had data are never touched - verified bit for bit - so
#' filling can only add information, never alter what was measured.
#'
#' Unlike the metrics here this is **not tiled**: GDAL fills the band as a
#' whole. On a raster far larger than memory, fill the tiles individually
#' before mosaicking.
#'
#' @section Tuning:
#' * `max_distance` is how far the interpolation will reach to find valid
#'   ground, in **map units**. Holes wider than twice this are left partly
#'   unfilled, which is usually what you want: a 2 km gap in coverage should
#'   not be invented from its edges.
#' * `smooth_iterations` runs 3x3 averaging passes over the filled cells,
#'   softening the spokes that inverse-distance interpolation can leave in a
#'   large hole. 0 is fine for scattered single-cell gaps; 2-3 helps on
#'   anything larger.
#'
#' @inheritParams rvt_svf
#' @param max_distance how far to search for valid cells, in **map units**
#'   (metres, normally; default 100)
#' @param smooth_iterations 3x3 smoothing passes applied to the filled cells
#'   only (default 0)
#' @return `out_path`, invisibly
#' @seealso [rvt_smooth()], which removes noise rather than gaps - run this
#'   first, since smoothing cannot sensibly cross a hole.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_fill(dem)
#' @export
rvt_fill <- function(dem, out_path = fs::file_temp(ext = "tif"),
                      max_distance = 100, smooth_iterations = 0,
                      threads = rvt_threads(), overwrite = FALSE) {
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  max_px <- .cells(max_distance, info$xres, "max_distance", "outer")
  smooth_iterations <- as.integer(smooth_iterations)
  if (is.na(smooth_iterations) || smooth_iterations < 0)
    stop("`smooth_iterations` must be zero or more", call. = FALSE)

  # fillNodata() edits a band in place, so it works on a scratch copy that is
  # then converted to the Cloud-Optimized form everything else here writes.
  scratch <- fs::file_temp(ext = "tif")
  on.exit(.rm_path(scratch), add = TRUE)
  gdalraster::translate(dem, scratch,
    cl_arg = c("-co", "TILED=YES", "-co", "COMPRESS=DEFLATE"), quiet = TRUE)
  gdalraster::fillNodata(scratch, band = 1L, max_dist = max_px,
                          smooth_iterations = smooth_iterations, quiet = TRUE)
  .finalize_cog(scratch, out_path, threads)
  invisible(out_path)
}

.smooth_tile <- function(tile, xres, yres, radius, norm_diff, iterations,
                          max_diff, pad, threads) {
  smooth_kernel(.pad_edge(tile, pad), pad, nrow(tile), ncol(tile), xres, yres,
                 radius, norm_diff, iterations, max_diff, threads)
}

#' Feature-preserving smoothing
#'
#' Removes noise from a surface while leaving breaks of slope sharp. An
#' ordinary averaging or Gaussian filter cannot tell a noisy cell from a
#' terrace edge and rounds off both, which is self-defeating here: the edges
#' are what most of these metrics exist to find.
#'
#' Worth reaching for when a DEM is speckled - gridded finer than the point
#' density really supports, or affected by flight-strip misalignment - and
#' before [rvt_curvature()], [rvt_log()] or [rvt_geomorphons()], all of which
#' amplify high-frequency noise by construction.
#'
#' @section How it works:
#' The algorithm is Sun, Rosin, Martin and Langbein's (2007), as used by
#' WhiteboxTools' `FeaturePreservingSmoothing`. It smooths the surface
#' *normals* rather than the elevations:
#'
#' * A unit normal is computed per cell, from the same central differences
#'   [rvt_slope()] uses.
#' * Each normal is replaced by a weighted average of those within `radius`,
#'   with the weight falling to zero once a neighbour's normal differs by more
#'   than `norm_diff` degrees. Across a break of slope the two sides do not
#'   see each other at all, so the edge survives, while noise - small random
#'   deviations well inside the threshold - averages away.
#' * Elevations are then nudged repeatedly towards agreement with those
#'   smoothed normals: each neighbour's plane implies a height for this cell,
#'   and the cell moves to their mean.
#'
#' That last step works on immediate neighbours only, not the whole `radius`
#' window - the radius belongs to the normal filter. Letting it reach further
#' lets cells on opposite sides of a scarp pull on each other and smears the
#' very edge the normal filter preserved.
#'
#' Total movement is capped at `max_diff`, so no cell can drift far from
#' measured ground however many iterations are run. Tiles overlap by the full
#' reach of the three stages, so the result does not depend on `tile_size`.
#'
#' @section Tuning:
#' * `norm_diff` is the one that matters. It is the angle, in degrees, beyond
#'   which two pieces of ground are treated as belonging to different
#'   surfaces. Low values (5-10) preserve almost every edge and remove only
#'   the finest noise; high values (25-40) smooth much harder and start
#'   rounding real breaks. The default 15 follows WhiteboxTools.
#' * `radius` sets how far the averaging reaches, in **map units**. Larger
#'   removes more noise and costs proportionally more - the window is
#'   circular, so cost grows with the square.
#' * `iterations` repeats only the elevation step, against normals that are
#'   filtered once. Three is usually enough; more mostly runs into `max_diff`.
#' * `max_diff` is a safety rail in elevation units. Lower it if you want to
#'   guarantee the surface stays close to what was measured.
#'
#' @inheritParams rvt_svf
#' @param radius how far the normal filter reaches, in **map units** (metres,
#'   normally; default 5)
#' @param norm_diff angle in degrees beyond which neighbouring ground is
#'   treated as a different surface and not averaged across (default 15)
#' @param iterations elevation-update passes (default 3)
#' @param max_diff largest change any cell may undergo, in elevation units
#'   (default 0.5)
#' @return `out_path`, invisibly
#' @references
#' Sun, X., Rosin, P.L., Martin, R.R. and Langbein, F.C. (2007) Fast and
#' effective feature-preserving mesh denoising. *IEEE Transactions on
#' Visualization and Computer Graphics* 13(5), 925-938.
#' \doi{10.1109/TVCG.2007.1065}
#' @seealso [rvt_fill()], which closes gaps rather than removing noise. Run it
#'   first: smoothing cannot sensibly average across a hole.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#'
#' # the usual preparation, as a pipeline
#' dem |> rvt_fill() |> rvt_smooth() |> rvt_curvature() |> rvt_plot()
#' @export
rvt_smooth <- function(dem, out_path = fs::file_temp(ext = "tif"),
                        radius = 5, norm_diff = 15, iterations = 3,
                        max_diff = 0.5,
                        tile_size = NULL, threads = rvt_threads(),
                        overwrite = FALSE, progress = FALSE) {
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  if (norm_diff <= 0 || norm_diff >= 90)
    stop("`norm_diff` must be between 0 and 90 degrees", call. = FALSE)
  iterations <- as.integer(iterations)
  if (is.na(iterations) || iterations < 0)
    stop("`iterations` must be zero or more", call. = FALSE)
  if (!is.finite(max_diff) || max_diff <= 0)
    stop("`max_diff` must be a positive elevation difference", call. = FALSE)
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  r_px <- .cells(radius, info$xres, "radius", "outer")
  # the reach of all three stages: normals (1 cell), the normal filter
  # (`radius`), then one cell per elevation pass, which works on immediate
  # neighbours only
  pad <- as.integer(r_px + iterations + 2L)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, pad)

  .process_tiled(dem, list(smooth = out_path), pad, tile_size,
                  function(tile, xres, yres)
                    list(smooth = .smooth_tile(tile, xres, yres, r_px,
                                                norm_diff, iterations, max_diff,
                                                pad, threads)),
                  progress = progress, threads = threads)
  invisible(out_path)
}
