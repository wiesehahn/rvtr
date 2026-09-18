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

## Every distance a user gives is in map units (metres, normally) and has to
## become a whole number of cells. Three rules, by what the distance is for:
##
##   "outer"  the far edge of a search window. Below one cell the raster simply
##            cannot answer the question, so this is an error rather than a
##            silent zero - which would otherwise produce a plausible-looking
##            result that means nothing (a sky-view factor of 1 everywhere).
##   "inner"  an inner radius or skip, where 0 is meaningful ("don't skip").
##   "step"   spacing between samples. Finer than the grid just means "every
##            cell", so it is clamped rather than refused.
##
## The halfway case is spelled out rather than left to round(), which rounds
## halves to *even*: round(2.5) is 2 while round(3.5) is 4, so a 5 m radius on
## a 2 m DEM and a 7 m radius on the same DEM would round in opposite
## directions.
.cells <- function(d, res, what, kind = c("outer", "inner", "step")) {
  kind <- match.arg(kind)
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d < 0)
    stop("`", what, "` must be a single non-negative distance in map units",
          call. = FALSE)

  # Tested against the request, not the rounded count: half a cell would round
  # *up* to one, which would answer a question the raster cannot resolve.
  if (kind == "outer" && d < res * (1 - 1e-9))
    stop(sprintf("`%s` = %s is smaller than one cell (%s); the finest this raster supports is %s.",
                  what, format(d), format(res), format(res)), call. = FALSE)

  px <- floor(d / res + 0.5)
  if (kind == "step") px <- max(1, px)
  px <- as.integer(px)

  # Say so only when the grid could not honour the request closely - silent for
  # the usual case of a round distance on a matching resolution.
  eff <- px * res
  if (d > 0 && abs(eff - d) > 0.1 * d)
    message(sprintf("`%s`: %s snapped to %d cell%s = %s (cell size %s)",
                     what, format(d), px, if (px == 1L) "" else "s",
                     format(eff), format(res)))
  px
}

## Every path a user hands in goes through here. fs normalises separators, so
## a Windows path written with backslashes reaches GDAL as something it can
## open, and expands "~", which GDAL does not do at all - `rvt_svf(dem,
## "~/svf.tif")` would otherwise create a directory called "~".
.as_path <- function(p) fs::path_expand(fs::path(p))

## Delete if present. fs::file_delete() errors on a missing file, which is the
## wrong behaviour in an on.exit() handler that may run after a failure.
.rm_path <- function(p) {
  p <- p[fs::file_exists(p)]
  if (length(p)) fs::file_delete(p)
  invisible(NULL)
}

.dem_info <- function(path, band = 1L) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  res <- ds$res()
  list(nx = ds$getRasterXSize(), ny = ds$getRasterYSize(),
       xres = res[1], yres = res[2],
       nodata = ds$getNoDataValue(band))
}

## Read a window from an already-open dataset. The tiling loop opens each
## source, pyramid level and output once and keeps the handle: reopening per
## window costs ~8 ms, which is nothing beside four big tiles but 0.5 s over
## 64 small ones, per input and per output.
.read_window_ds <- function(ds, x0, y0, ncols, nrows, nodata, band = 1L) {
  vals <- ds$read(band, x0, y0, ncols, nrows, ncols, nrows)
  m <- matrix(as.numeric(vals), nrow = nrows, ncol = ncols, byrow = TRUE)
  if (!is.na(nodata)) m[m == nodata] <- NA
  m
}

.read_window <- function(path, x0, y0, ncols, nrows, nodata, band = 1L) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  .read_window_ds(ds, x0, y0, ncols, nrows, nodata, band)
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
                                xres, yres, oversample = 3, angles = NULL) {
  if (is.null(angles))
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

## Offsets indexed by *compass bearing* rather than by the counterclockwise-
## from-east convention above: direction d points at bearing d * 360/n. The
## daylight kernel needs this so a sun azimuth maps onto a direction index by
## arithmetic instead of a lookup.
.azimuth_offsets <- function(num_directions, radius_max, xres, yres,
                              oversample = 3, rmin = 1) {
  bearings <- (360 / num_directions) * (0:(num_directions - 1))
  .direction_offsets(num_directions, radius_max, rmin, xres, yres, oversample,
                      angles = ((90 - bearings) %% 360) * pi / 180)
}

## Minimum search radius grows with the noise-removal level, matching RVT's
## svf_noise (0 / 10 / 20 / 40 percent of radius_max).
.radius_min <- function(radius_max, noise_removal) {
  max(1, round(radius_max * c(0, 0.10, 0.20, 0.40)[noise_removal + 1]))
}

## Plain tiled GeoTIFF used as scratch space while tiles are written one
## window at a time (see .process_tiled()) - not the final output, see
## .finalize_cog().
## `dtype` is Float32 for data and Byte for display images, which are already
## 0-255 by the time they arrive: a Float32 scratch there pushes four times the
## bytes through disk for nothing.
##
## **Uncompressed on purpose.** This file is written once, read once by the
## finalizer and deleted. Compressing it cost more than it saved: measured on
## 4000x4000 Float32, DEFLATE+predictor at ZLEVEL 9 took 0.24 s to create and
## 1.10 s to write against 0.05 s and 0.33 s uncompressed, and the COG step
## then read it slightly faster too (0.65 s against 0.71 s) - about 1 s of the
## ~3.3 s a cheap metric costs end to end. ZLEVEL 1 was no faster than 9, so
## the cost is the predictor plus deflate itself, not the level. The trade is
## a larger temporary file, 4 bytes per cell rather than ~3.
.create_scratch <- function(src, path, tile_size, nx, ny, nbands = 1L,
                             nodata_value = -9999, dtype = "Float32") {
  # Blocks stay aligned to the write tiles: it no longer matters for
  # compression, but a window that straddles blocks still costs extra I/O.
  bs <- max(16, (min(tile_size, nx, ny) %/% 16) * 16)
  gdalraster::rasterFromRaster(
    srcfile = src, dstfile = path, nbands = nbands, dtName = dtype,
    init = if (is.null(nodata_value)) 0 else nodata_value,
    options = c("COMPRESS=NONE", "TILED=YES",
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
## Overviews: MODE for categorical results (geomorphons class codes), CUBIC
## otherwise. Averaging classes invents codes that aren't classes - measured on
## geomorphons, 36 distinct non-integer overview values against 9 whole classes
## with MODE. CUBIC was checked beside NoData holes before adopting it: same
## NoData count as AVERAGE and no values outside the data's range.
## OVERVIEWS=IGNORE_EXISTING because the COG writer otherwise *copies* a
## source's existing overviews and ignores RESAMPLING entirely (the scratch
## file has none, but the rule is cheap to keep uniform). PREDICTOR=FLOATING_POINT
## is the COG driver's name for what the plain GTiff driver above calls
## PREDICTOR=3 - same TIFF predictor tag, different driver, different string.
.finalize_cog <- function(scratch_path, out_path, threads, categorical = FALSE) {
  gdalraster::translate(scratch_path, out_path,
    cl_arg = c("-of", "COG",
                "-co", "COMPRESS=DEFLATE",
                "-co", "PREDICTOR=FLOATING_POINT",
                "-co", paste0("RESAMPLING=", if (categorical) "MODE" else "CUBIC"),
                "-co", "OVERVIEWS=IGNORE_EXISTING",
                "-co", paste0("NUM_THREADS=", threads)),
    quiet = TRUE)
  .rm_path(scratch_path)
  invisible(out_path)
}

## `bands` is a list of equally sized matrices, one per output band. The
## dataset is opened by the caller and stays open across tiles; NoData is set
## once when it is opened, not on every window.
.write_window_ds <- function(ds, bands, x0, y0, nodata_value = -9999) {
  fill <- if (is.null(nodata_value)) 0 else nodata_value
  for (b in seq_along(bands)) {
    mat <- bands[[b]]
    mat[is.na(mat)] <- fill
    ds$write(b, x0, y0, ncol(mat), nrow(mat), as.vector(t(mat)))
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
##
## `aux` adds further, coarser rasters read alongside each tile - the pyramid
## levels of .pyramid_plan(). Each entry is `list(path, fac, overlap)`, where
## `fac` is how many full-resolution cells one of its cells covers. For every
## tile the matching window is read from each, widened by that level's own
## overlap, and `fun` is called as `fun(tile, xres, yres, ctx)` with `ctx`
## carrying those windows plus the tile's global pixel offset (which the
## kernel needs to map an output cell onto a coarse cell).
##
## The aux windows are deliberately *not* padded. A window always spans
## [first coarse cell of the tile - overlap, last + overlap], so a sample can
## only fall outside it where the raster itself ends - which means clamping to
## the window is clamping to the raster, and stays identical however the
## raster is tiled.
.process_tiled <- function(src, out_paths, overlap, tile_size, fun,
                            band = 1L, progress = FALSE, threads = rvt_threads(),
                            nbands = NULL, aux = NULL, categorical = FALSE,
                            image = FALSE, quality = 90) {
  info <- .dem_info(src, band)
  nx <- info$nx; ny <- info$ny

  # An output ending in .webp/.jpg is a picture, finished by .finalize_image().
  # Only callers that produce 0-255 display bands may ask for one; a metric's
  # Float32 values cast to Byte would be silent garbage. Size limits are
  # checked here, before a single tile is computed.
  formats <- lapply(out_paths, .image_format)
  for (nm in names(out_paths)) {
    f <- formats[[nm]]
    if (is.null(f)) next
    if (!isTRUE(image))
      stop(paste("This function writes a GeoTIFF (`.tif`). For a picture, pipe",
                 "its result into rvt_image()."), call. = FALSE)
    .image_check_size(nx, ny, f)
  }

  if (!is.null(aux)) aux <- lapply(aux, function(a) {
    ai <- .dem_info(a$path, if (is.null(a$band)) 1L else a$band)
    c(a, list(nx = ai$nx, ny = ai$ny, nodata = ai$nodata))
  })

  nb <- stats::setNames(rep(1L, length(out_paths)), names(out_paths))
  if (!is.null(nbands)) nb[names(nbands)] <- as.integer(nbands)

  # A display image is already 0-255 when it arrives, so its scratch is Byte
  # and carries no NoData: holes are handled by the alpha band (WebP) or come
  # out black (JPEG).
  scratch_paths <- lapply(out_paths, function(p) fs::file_temp(ext = "tif"))
  nodata_values <- lapply(formats, function(f) if (is.null(f)) -9999 else NULL)
  for (nm in names(out_paths))
    .create_scratch(src, scratch_paths[[nm]], tile_size, nx, ny, nb[[nm]],
                    nodata_value = nodata_values[[nm]],
                    dtype = if (is.null(formats[[nm]])) "Float32" else "Byte")

  # One open handle per source, pyramid level and output for the whole loop
  src_ds <- methods::new(gdalraster::GDALRaster, src, read_only = TRUE)
  aux_ds <- lapply(aux, function(a)
    methods::new(gdalraster::GDALRaster, a$path, read_only = TRUE))
  out_ds <- lapply(scratch_paths, function(p)
    methods::new(gdalraster::GDALRaster, p, read_only = FALSE))
  on.exit({
    try(src_ds$close(), silent = TRUE)
    for (d in aux_ds) try(d$close(), silent = TRUE)
    for (d in out_ds) try(d$close(), silent = TRUE)
  }, add = TRUE)
  for (nm in names(out_paths)) {
    nd <- nodata_values[[nm]]
    if (!is.null(nd)) for (b in seq_len(nb[[nm]])) out_ds[[nm]]$setNoDataValue(b, nd)
  }

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

      tile <- .read_window_ds(src_ds, x0 - left, y0 - top,
                               cols + left + right, rows + top + bottom,
                               info$nodata, band)

      res <- if (is.null(aux)) {
        fun(tile, info$xres, info$yres)
      } else {
        # the kernels compute every cell of the tile, overlap band included,
        # so the coarse window has to span the tile's full extent - not just
        # the part that gets written
        gx0 <- x0 - left; gy0 <- y0 - top
        gx1 <- gx0 + (cols + left + right) - 1L
        gy1 <- gy0 + (rows + top + bottom) - 1L
        wins <- lapply(seq_along(aux), function(i) {
          a <- aux[[i]]
          cx0 <- max(0L, (gx0 %/% a$fac) - a$overlap)
          cy0 <- max(0L, (gy0 %/% a$fac) - a$overlap)
          cx1 <- min(a$nx - 1L, (gx1 %/% a$fac) + a$overlap)
          cy1 <- min(a$ny - 1L, (gy1 %/% a$fac) + a$overlap)
          list(m = .read_window_ds(aux_ds[[i]], cx0, cy0,
                                    cx1 - cx0 + 1L, cy1 - cy0 + 1L, a$nodata,
                                    if (is.null(a$band)) 1L else a$band),
                cx0 = as.integer(cx0), cy0 = as.integer(cy0),
                fac = as.integer(a$fac))
        })
        fun(tile, info$xres, info$yres,
            list(aux = wins, x0 = x0 - left, y0 = y0 - top))
      }

      for (nm in names(out_paths)) {
        m <- res[[nm]]
        if (!is.list(m)) m <- list(m)
        .write_window_ds(out_ds[[nm]],
                          lapply(m, function(b)
                            b[(top + 1):(top + rows), (left + 1):(left + cols), drop = FALSE]),
                          x0, y0, nodata_values[[nm]])
      }

      if (progress) message(sprintf("  tile %d/%d", k, n_tiles))
    }
  }

  # close every handle before the finalizers read the scratch files back
  src_ds$close()
  for (d in aux_ds) d$close()
  for (d in out_ds) d$close()

  for (nm in names(out_paths)) {
    if (is.null(formats[[nm]]))
      .finalize_cog(scratch_paths[[nm]], out_paths[[nm]], threads,
                    categorical = categorical)
    else
      .finalize_image(scratch_paths[[nm]], out_paths[[nm]], formats[[nm]], quality)
  }

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
