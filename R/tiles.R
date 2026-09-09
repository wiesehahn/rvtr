## Internal machinery: raster I/O, tiling, and the horizon-search geometry.
## Everything here works on a single GDAL dataset path - which may be a plain
## GeoTIFF or a .vrt standing in for a whole mosaic (see rvt_mosaic()), so the
## tiling engine never needs to know whether it is looking at one file or a
## thousand.

#' Number of threads used by the C++ kernels
#'
#' Defaults to all detected cores. Set with `rvt_threads(n)` or by setting
#' `options(rvtr.threads = n)`.
#'
#' @param n number of threads, or NULL to just query the current setting
#' @return the number of threads in effect
#' @examples
#' rvt_threads()
#' @export
rvt_threads <- function(n = NULL) {
  if (!is.null(n)) options(rvtr.threads = as.integer(n))
  getOption("rvtr.threads", default = max(1L, parallel::detectCores()))
}

.dem_info <- function(path, band = 1L) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  res <- ds$res()
  list(nx = ds$getRasterXSize(), ny = ds$getRasterYSize(),
       xres = res[1], yres = res[2],
       nodata = ds$getNoDataValue(band))
}

.read_window <- function(path, x0, y0, ncols, nrows, nodata, band = 1L) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  vals <- ds$read(band, x0, y0, ncols, nrows, ncols, nrows)
  m <- matrix(as.numeric(vals), nrow = nrows, ncol = ncols, byrow = TRUE)
  if (!is.na(nodata)) m[m == nodata] <- NA
  m
}

## Replicate the edge cell outwards (numpy's mode="edge"), the standard edge
## handling for finite-difference derivatives. Used for slope/aspect.
.pad_edge <- function(mat, pad) {
  nr <- nrow(mat); nc <- ncol(mat)
  mat <- rbind(mat[rep(1L, pad), , drop = FALSE], mat,
               mat[rep(nr, pad), , drop = FALSE])
  cbind(mat[, rep(1L, pad), drop = FALSE], mat,
        mat[, rep(nc, pad), drop = FALSE])
}

## Mirror-pad *repeating* the edge cell (numpy's mode="symmetric"), and unlike
## .pad_reflect() it copes with `pad` larger than the matrix by folding
## repeatedly. MSTP needs this: its windows reach 2023 px, far beyond a small
## raster, and edge-replicate padding there would fill the window with copies
## of one value, collapsing the standard deviation MSTP divides by.
.mirror_index <- function(n, pad) {
  i <- seq_len(n + 2 * pad) - pad - 1L      # 0-based source index, may be <0
  period <- 2L * n
  i <- i %% period
  ifelse(i >= n, period - 1L - i, i) + 1L   # fold, back to 1-based
}

.pad_symmetric <- function(mat, pad) {
  mat[.mirror_index(nrow(mat), pad), .mirror_index(ncol(mat), pad), drop = FALSE]
}

## Mirror-pad without repeating the edge cell (numpy's mode="reflect"), used
## to extend the raster beyond its own edge for the horizon search - it
## extrapolates terrain more plausibly there than a flat edge-replicate would.
.pad_reflect <- function(mat, pad) {
  nr <- nrow(mat); nc <- ncol(mat)
  if (pad >= nr || pad >= nc)
    stop("tile smaller than the search radius; increase tile_size", call. = FALSE)
  mat <- rbind(mat[(pad + 1):2, , drop = FALSE],
               mat,
               mat[(nr - 1):(nr - pad), , drop = FALSE])
  cbind(mat[, (pad + 1):2, drop = FALSE],
        mat,
        mat[, (nc - 1):(nc - pad), drop = FALSE])
}

## Integer pixel offsets to search per direction, deduplicated and sorted by
## ascending distance, flattened into the layout the C++ kernel expects.
## Angles run counterclockwise from east; row shifts are negated because the
## raster row index increases southward.
.direction_offsets <- function(num_directions, radius_max, radius_min,
                                xres, yres, oversample = 3) {
  angles <- (2 * pi / num_directions) * (0:(num_directions - 1))
  radii <- seq(radius_min, radius_max, by = 1 / oversample)

  per_dir <- lapply(angles, function(a) {
    dx <- round(cos(a) * radii)
    dy <- round(-sin(a) * radii)
    keep <- !duplicated(paste(dx, dy)) & !(dx == 0 & dy == 0)
    dx <- dx[keep]; dy <- dy[keep]
    ord <- order(sqrt(dx^2 + dy^2))
    list(dx = dx[ord], dy = dy[ord])
  })

  dx <- unlist(lapply(per_dir, `[[`, "dx"))
  dy <- unlist(lapply(per_dir, `[[`, "dy"))
  n_per <- vapply(per_dir, function(o) length(o$dx), integer(1))
  ends <- cumsum(n_per)
  list(dx = as.integer(dx), dy = as.integer(dy),
       dist = sqrt((dx * xres)^2 + (dy * yres)^2),   # ground distance
       starts = as.integer(c(0L, utils::head(ends, -1))),
       ends = as.integer(ends),
       angles = angles)
}

## Minimum search radius grows with the noise-removal level, matching RVT's
## svf_noise (0 / 10 / 20 / 40 percent of radius_max).
.radius_min <- function(radius_max, noise_removal) {
  max(1, round(radius_max * c(0, 0.10, 0.20, 0.40)[noise_removal + 1]))
}

## Plain tiled GeoTIFF used as scratch space while tiles are written one
## window at a time (see .process_tiled()) - not the final output, see
## .finalize_cog().
.create_scratch <- function(src, path, tile_size, nx, ny, nbands = 1L,
                             nodata_value = -9999) {
  # Internal block size is aligned to the write tiles: if blocks straddle tile
  # boundaries GDAL decompresses and recompresses the same block once per
  # touching tile, which inflates the file several-fold.
  bs <- max(16, (min(tile_size, nx, ny) %/% 16) * 16)
  gdalraster::rasterFromRaster(
    srcfile = src, dstfile = path, nbands = nbands, dtName = "Float32",
    init = nodata_value,
    options = c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=9", "TILED=YES",
                 paste0("BLOCKXSIZE=", bs), paste0("BLOCKYSIZE=", bs)),
    quiet = TRUE)
  invisible(path)
}

## Convert the finished scratch GeoTIFF into a Cloud-Optimized GeoTIFF at the
## path the caller actually asked for. GDAL's COG driver only supports
## CreateCopy() (it needs a complete source to compute overviews from), not
## incremental RasterIO writes - hence writing to a plain tiled GeoTIFF above
## and converting here, rather than writing the COG directly.
##
## AVERAGE resampling for overviews: it's NoData-aware (pixels flagged via
## setNoDataValue() are excluded from the mean, not averaged in as if they
## were real low values) and appropriate for continuous data like SVF/openness
## /VAT - nearest-neighbour (the driver's default) would look blocky in a
## quick low-res preview. PREDICTOR=FLOATING_POINT is the COG driver's name
## for what the plain GTiff driver above calls PREDICTOR=3 - same TIFF
## predictor tag, different driver, different string.
.finalize_cog <- function(scratch_path, out_path, threads) {
  gdalraster::translate(scratch_path, out_path,
    cl_arg = c("-of", "COG",
                "-co", "COMPRESS=DEFLATE",
                "-co", "PREDICTOR=FLOATING_POINT",
                "-co", "RESAMPLING=AVERAGE",
                "-co", paste0("NUM_THREADS=", threads)),
    quiet = TRUE)
  unlink(scratch_path)
  invisible(out_path)
}

## `bands` is a list of equally sized matrices, one per output band.
.write_window <- function(path, bands, x0, y0, nodata_value = -9999) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = FALSE)
  on.exit(ds$close())
  for (b in seq_along(bands)) {
    mat <- bands[[b]]
    mat[is.na(mat)] <- nodata_value
    ds$write(b, x0, y0, ncol(mat), nrow(mat), as.vector(t(mat)))
    ds$setNoDataValue(b, nodata_value)
  }
  invisible(NULL)
}

## Tile-by-tile driver. `fun(tile, xres, yres)` receives a tile *including* its
## overlap border and must return a named list, one entry per output: either a
## matrix, or a list of matrices for a multi-band output. The overlap is
## cropped here before writing. Only one tile is held in memory at a time, so
## input size is bounded by disk, not RAM.
##
## `nbands` says how many bands each output has, and has to be known up front
## because the scratch file is created before the first tile runs. Named to
## match `out_paths`; anything not listed is single-band.
##
## Outputs are written to scratch files during tiling, then converted to
## Cloud-Optimized GeoTIFFs at the requested paths once complete (see
## .finalize_cog()) - so a caller can cheaply read a low-res overview of a
## huge result (rvt_plot() does exactly this) or a full-resolution window of
## it, without touching the rest of the file.
.process_tiled <- function(src, out_paths, overlap, tile_size, fun,
                            band = 1L, progress = FALSE, threads = rvt_threads(),
                            nbands = NULL) {
  info <- .dem_info(src, band)
  nx <- info$nx; ny <- info$ny

  nb <- stats::setNames(rep(1L, length(out_paths)), names(out_paths))
  if (!is.null(nbands)) nb[names(nbands)] <- as.integer(nbands)

  scratch_paths <- lapply(out_paths, function(p) tempfile(fileext = ".tif"))
  for (nm in names(out_paths))
    .create_scratch(src, scratch_paths[[nm]], tile_size, nx, ny, nb[[nm]])

  x0s <- seq(0, nx - 1, by = tile_size)
  y0s <- seq(0, ny - 1, by = tile_size)
  n_tiles <- length(x0s) * length(y0s)
  k <- 0L

  for (y0 in y0s) {
    for (x0 in x0s) {
      k <- k + 1L
      cols <- min(tile_size, nx - x0)
      rows <- min(tile_size, ny - y0)

      left   <- min(overlap, x0)
      top    <- min(overlap, y0)
      right  <- min(overlap, nx - x0 - cols)
      bottom <- min(overlap, ny - y0 - rows)

      tile <- .read_window(src, x0 - left, y0 - top,
                            cols + left + right, rows + top + bottom,
                            info$nodata, band)

      res <- fun(tile, info$xres, info$yres)

      for (nm in names(out_paths)) {
        m <- res[[nm]]
        if (!is.list(m)) m <- list(m)
        .write_window(scratch_paths[[nm]],
                       lapply(m, function(b)
                         b[(top + 1):(top + rows), (left + 1):(left + cols), drop = FALSE]),
                       x0, y0)
      }

      if (progress) message(sprintf("  tile %d/%d", k, n_tiles))
    }
  }

  for (nm in names(out_paths)) .finalize_cog(scratch_paths[[nm]], out_paths[[nm]], threads)

  invisible(out_paths)
}

## Choose a tile size that keeps peak memory near `target_mb` per matrix while
## staying comfortably larger than the overlap band.
.auto_tile_size <- function(nx, ny, overlap, target_mb = 256) {
  cells <- target_mb * 1024^2 / 8
  ts <- floor(sqrt(cells))
  ts <- max(ts, 4L * overlap + 16L)
  min(ts, max(nx, ny))
}
