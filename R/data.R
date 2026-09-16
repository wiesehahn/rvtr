## Helpers shared by the rvt_data_*() fetchers.

## Snap an extent outward to whole multiples of `res`, so rasters fetched for
## the same place at the same resolution share a grid.
.data_snap <- function(extent, res) {
  r <- function(v) round(v, 6)
  c(r(floor(extent[1] / res + 1e-9) * res), r(floor(extent[2] / res + 1e-9) * res),
    r(ceiling(extent[3] / res - 1e-9) * res), r(ceiling(extent[4] / res - 1e-9) * res))
}

## Copy a (usually virtual, remote) raster into a Cloud-Optimized GeoTIFF.
## Written to a hidden file beside the target, then moved into place, so an
## interrupted download never leaves a truncated file that looks valid.
## `metadata` is a named character vector. Returns FALSE if GDAL failed, so the
## caller can say what it was fetching.
.data_write_cog <- function(src, out_path, predictor, metadata, threads) {
  fs::dir_create(fs::path_dir(out_path))
  tmp <- fs::path(fs::path_dir(out_path), paste0(".", fs::path_file(out_path), ".part"))
  .rm_path(tmp)
  on.exit(.rm_path(tmp), add = TRUE)
  mo <- as.vector(rbind("-mo", paste0(names(metadata), "=", metadata)))
  ok <- tryCatch({
    gdalraster::translate(src, tmp,
      cl_arg = c("-of", "COG",
                 "-co", "COMPRESS=DEFLATE",
                 "-co", paste0("PREDICTOR=", predictor),
                 "-co", "RESAMPLING=CUBIC",
                 "-co", "OVERVIEWS=IGNORE_EXISTING",
                 "-co", paste0("NUM_THREADS=", threads),
                 mo),
      quiet = TRUE)
    fs::file_exists(tmp)
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) return(FALSE)
  .rm_path(out_path)
  fs::file_move(tmp, out_path)
  TRUE
}
