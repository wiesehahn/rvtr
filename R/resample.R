## Changing the scale of analysis by coarsening the grid.

## Which resampling algorithms each GDAL utility can actually perform.
## gdal_translate goes through RasterIO, which knows only the first set and
## silently substitutes nearest neighbour for anything else; gdalwarp knows
## both. Keep these in step with GDAL if the lists ever grow.
.rasterio_methods <- c("near", "bilinear", "cubic", "cubicspline", "lanczos",
                        "average", "rms", "mode", "gauss")
.warp_only_methods <- c("max", "min", "med", "q1", "q3", "sum")

#' Resample a raster to a coarser resolution
#'
#' Nearly every metric in this package answers a different question at a
#' different resolution: on a 0.25 m DTM [rvt_curvature()] picks out cart ruts
#' and tree-throw pits, on a 2 m version of the same ground it picks out
#' terrace edges and gully walls. This coarsens a raster so you can ask at the
#' scale you actually mean.
#'
#' It only ever goes **coarser**. Detail that was never measured cannot be
#' invented, so asking for a finer resolution than the source is an error
#' rather than a silently smooth upsampling.
#'
#' @section How it works:
#' A wrapper around GDAL's `gdal_translate -tr` (via
#' [gdalraster::translate()]) that writes a real Cloud-Optimized GeoTIFF
#' rather than a lazy VRT reference. Writing the file is the point:
#'
#' * The values get frozen once. GDAL's `average` resampling at a non-integer
#'   scale factor gives *different elevations depending on how the raster is
#'   read* - reading the bundled tile at 2.5 m in one go versus in 137-pixel
#'   windows disagrees by up to 0.53 m. Since every metric here walks its
#'   input in tiles, a lazy reference would make results depend on
#'   `tile_size`, which is the one thing this package guarantees never
#'   happens.
#' * Reads afterwards are about 10x faster than re-resampling on every read,
#'   and you will normally run several metrics on the same coarsened DEM.
#'
#' Overviews on the source are explicitly ignored (`-ovr NONE`). GDAL would
#' otherwise quietly satisfy a coarsening request from whatever overviews the
#' file happens to carry, which silently overrides `method` - asking for
#' `"mode"` on a raster with averaged overviews returns averages, not classes.
#'
#' @section Choosing a method:
#' `"cubic"` (the default) suits continuous surfaces - elevation, and any of
#' the continuous products here. It keeps the shape of the terrain better than
#' a plain block average, and it is **NoData-safe**: a cell whose neighbourhood
#' is incomplete becomes NoData rather than being blended with the fill value.
#' The price is that holes grow a little - on the bundled tile, coarsening 1 m
#' to 4 m leaves 141 NoData cells against `"average"`'s 9.
#'
#' `"mode"` (majority) is the one to use for **categorical** rasters - the
#' class codes from [rvt_geomorphons()]. Each coarse cell takes the commonest
#' class beneath it, so the result still contains only whole class codes.
#' `"near"` also preserves codes but picks one arbitrary cell per block rather
#' than the commonest, so prefer `"mode"` unless you specifically want
#' subsampling.
#'
#' **Avoid `"average"` on anything this package wrote.** GDAL's average does
#' not exclude the NoData value from the mean, and every output here uses
#' -9999 - so cells beside a hole come out around -9000 instead of NoData.
#' It is safe only on rasters whose NoData is NaN (which the bundled DEM
#' happens to use). A warning is issued if you ask for it on a raster with a
#' numeric NoData value.
#'
#' `"max"` is what the multi-resolution horizon search ([rvt_reach]) builds its
#' coarse levels with, because it is the only choice that keeps a distant
#' skyline intact: averaging flattens a narrow obstruction and nearest
#' neighbour can miss it altogether.
#'
#' Any other algorithm GDAL accepts is passed through (`"bilinear"`,
#' `"lanczos"`, `"min"`, `"med"`, `"rms"`, `"q1"`, `"q3"`, `"sum"`, ...) if you
#' know you want it. Note `"max"`, `"min"`, `"med"`, `"q1"`, `"q3"` and
#' `"sum"` are done with `gdalwarp` rather than `gdal_translate`, which cannot
#' perform them, and so need the raster to carry a coordinate reference
#' system.
#'
#' @section Which scale knob to reach for:
#' Resolution is only half of scale - the other half is how far each metric
#' looks, and the two are not interchangeable. Coarsening removes the
#' microrelief *before* the metric sees it; widening a search radius keeps the
#' microrelief and measures over more ground.
#'
#' * [rvt_slope()], [rvt_aspect()], [rvt_curvature()] and [rvt_log()] work off
#'   a fixed 3x3 window and have no radius to widen - resampling is their only
#'   scale control.
#' * [rvt_slrm()], [rvt_msrm()], [rvt_dev()] and [rvt_mstp()] cost the same per
#'   pixel whatever their radius, so widening the radius is free and keeps more
#'   information. Reach for that first.
#' * [rvt_svf()], [rvt_openness()], [rvt_sky_illumination()] and [rvt_shadow()]
#'   cost roughly cells x radius, so coarsening pays twice over: fewer cells
#'   *and* shorter rays for the same reach in metres. Measured on a 4000x4000
#'   0.25 m tile, sky-view factor over a 10 m reach takes 15.8 s; on the 1 m
#'   version of the same ground, 0.6 s. That is often the difference between
#'   practical and not on a large survey.
#'
#' Every distance in this package is in **map units**, so coarsening does not
#' silently change what a parameter means: `rvt_svf(dem, reach = 10)` looks
#' 10 m out whatever the resolution. What changes is how finely that 10 m is
#' sampled - which is the point of coarsening.
#'
#' @param src path to a raster - a DEM, or any single- or multi-band result
#'   from this package. Several paths are mosaicked first, as elsewhere.
#' @param res target cell size in the raster's own map units (usually metres),
#'   either one number for square cells or two for `c(x, y)`
#' @param method resampling algorithm: `"cubic"` (default, continuous data),
#'   `"mode"` (categorical), or any other algorithm
#'   `gdal_translate -r` accepts. See **Choosing a method**.
#' @param out_path where to write the result (default: a temporary file)
#' @param threads threads for the output compression
#' @param overwrite overwrite `out_path` if it already exists
#' @return `out_path`, invisibly - so it pipes straight into any metric.
#' @seealso [rvt_mosaic()], which joins tiles without changing resolution.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#'
#' # ask a landscape-scale question of a 1 m DEM
#' dem |> rvt_resample(4) |> rvt_curvature() |> rvt_plot()
#' @export
rvt_resample <- function(src, res, method = "cubic",
                          out_path = fs::file_temp(ext = "tif"),
                          threads = rvt_threads(), overwrite = FALSE) {
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  src <- rvt_mosaic(src)

  res <- as.numeric(res)
  if (length(res) == 1L) res <- c(res, res)
  if (length(res) != 2L || anyNA(res) || any(res <= 0))
    stop("`res` must be one or two positive numbers", call. = FALSE)

  info <- .dem_info(src)
  # A relative tolerance, so asking for exactly the source resolution is a
  # no-op copy rather than a rounding-error failure.
  if (res[1] < info$xres * (1 - 1e-9) || res[2] < info$yres * (1 - 1e-9))
    stop(sprintf(paste("`res` (%s) is finer than the source (%s x %s).",
                        "Resampling can only coarsen - detail that was never",
                        "measured cannot be invented."),
                  paste(format(res), collapse = " x "),
                  format(info$xres), format(info$yres)), call. = FALSE)

  ds <- methods::new(gdalraster::GDALRaster, src, read_only = TRUE)
  dtype <- ds$getDataTypeName(1L)
  srs <- ds$getProjectionRef()
  ds$close()

  if (method == "average" && !is.na(info$nodata))
    warning(sprintf(paste("`method = \"average\"` does not exclude the NoData",
                           "value (%s) from the mean, so cells beside a hole",
                           "will be contaminated. Use \"cubic\" for continuous",
                           "data or \"mode\" for classes."),
                     format(info$nodata)), call. = FALSE)

  # Overviews are built with the same intent as the resampling itself:
  # averaging a classification would invent codes that aren't classes.
  ovr <- if (method %in% c("mode", "near")) "MODE" else "AVERAGE"
  predictor <- if (dtype %in% c("Float32", "Float64")) "FLOATING_POINT"
                else if (dtype == "Byte") "NONE" else "STANDARD"

  common <- c("-tr", format(res[1], scientific = FALSE),
                     format(res[2], scientific = FALSE),
              "-r", method,
              "-ovr", "NONE",
              "-of", "COG",
              "-co", "COMPRESS=DEFLATE",
              "-co", paste0("PREDICTOR=", predictor),
              "-co", paste0("RESAMPLING=", ovr),
              "-co", paste0("NUM_THREADS=", threads))

  if (method %in% .warp_only_methods) {
    # gdal_translate resamples through RasterIO, which supports only part of
    # the algorithm list and - this is the trap - **silently falls back to
    # nearest neighbour** for the rest, warning only at GDAL's level 6. Asking
    # translate for "max" gave output bit-identical to "near". gdalwarp
    # implements these properly, so they go that way instead.
    if (!nzchar(srs))
      stop(sprintf(paste("`method = \"%s\"` needs the raster to carry a",
                          "coordinate reference system (it is done with",
                          "gdalwarp, which requires one). Assign a CRS, or use",
                          "one of: %s."),
                    method, paste(.rasterio_methods, collapse = ", ")),
            call. = FALSE)
    gdalraster::warp(src, out_path, t_srs = srs, cl_arg = common, quiet = TRUE)
  } else {
    gdalraster::translate(src, out_path, cl_arg = common, quiet = TRUE)
  }

  invisible(out_path)
}
