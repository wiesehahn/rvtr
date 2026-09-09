## Sky-view factor and positive openness.

## Run the fused kernel over one tile (which already carries its overlap
## border). The reflect padding added here only affects the outermost band,
## which is either cropped away by the tiling engine or - at the true raster
## edge - is exactly the extension RVT uses.
.horizon_tile <- function(tile, xres, yres, num_directions, radius_max,
                           noise_removal, want_svf, want_opns, threads) {
  off <- .direction_offsets(num_directions, radius_max,
                             .radius_min(radius_max, noise_removal), xres, yres)
  pad <- as.integer(radius_max + 1L)
  padded <- .pad_reflect(tile, pad)
  out <- horizon_svf_opns(padded, pad, nrow(tile), ncol(tile),
                           off$dx, off$dy, off$dist, off$starts, off$ends,
                           want_svf, want_opns, threads)
  list(svf = out$svf, opns = out$opns)
}

.horizon_run <- function(dem, out_path, want_svf, want_opns,
                          num_directions, radius_max, noise_removal,
                          tile_size, threads, overwrite, progress) {
  if (!overwrite && file.exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  overlap <- as.integer(radius_max + 1L)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

  nm <- if (want_svf) "svf" else "opns"
  .process_tiled(dem, stats::setNames(list(out_path), nm), overlap, tile_size,
                  function(tile, xres, yres)
                    .horizon_tile(tile, xres, yres, num_directions, radius_max,
                                   noise_removal, want_svf, want_opns, threads),
                  progress = progress, threads = threads)
  invisible(out_path)
}

#' Sky-view factor
#'
#' Portion of the sky visible from each cell, in the range 0-1. Follows the
#' Relief Visualization Toolbox implementation (Zaksek et al. 2011).
#'
#' @param dem path to a DEM, or a vector of paths forming a mosaic (see
#'   [rvt_mosaic()]); mosaics are read across file boundaries, so tiles do not
#'   produce seams
#' @param out_path output GeoTIFF path (default: a temp file)
#' @param num_directions number of search directions (RVT default 16)
#' @param radius_max maximum search radius in *pixels* (RVT default 10)
#' @param noise_removal 0-3, raises the minimum search radius to ignore
#'   near-cell noise (RVT svf_noise: 0/10/20/40 percent of `radius_max`)
#' @param tile_size tile edge in pixels; NULL picks a size targeting roughly
#'   256 MB per working matrix
#' @param threads C++ threads (see [rvt_threads()])
#' @param overwrite recompute even if `out_path` exists
#' @param progress report per-tile progress
#' @return `out_path`, invisibly. `dem` is the first argument, so this is
#'   pipe-friendly: `dem |> rvt_svf("svf.tif")`.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_svf(dem)
#' @export
rvt_svf <- function(dem, out_path = tempfile(fileext = ".tif"),
                num_directions = 16, radius_max = 10, noise_removal = 0,
                tile_size = NULL, threads = rvt_threads(),
                overwrite = FALSE, progress = FALSE) {
  .horizon_run(dem, out_path, TRUE, FALSE, num_directions, radius_max,
                noise_removal, tile_size, threads, overwrite, progress)
}

#' Positive openness
#'
#' Mean zenith angle of the horizon, in degrees.
#'
#' @inheritParams rvt_svf
#' @return `out_path`, invisibly
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_openness(dem)
#' @export
rvt_openness <- function(dem, out_path = tempfile(fileext = ".tif"),
                     num_directions = 16, radius_max = 10, noise_removal = 0,
                     tile_size = NULL, threads = rvt_threads(),
                     overwrite = FALSE, progress = FALSE) {
  .horizon_run(dem, out_path, FALSE, TRUE, num_directions, radius_max,
                noise_removal, tile_size, threads, overwrite, progress)
}
