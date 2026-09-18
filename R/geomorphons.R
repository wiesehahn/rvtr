## Geomorphons landform classification.

#' Geomorphon landform classes
#'
#' The ten classes [rvt_geomorphons()] assigns, in code order: a named integer
#' vector mapping class name to the value stored in the raster.
#'
#' @format Named integer vector of length 10.
#' @examples
#' rvt_geomorphon_classes
#' names(rvt_geomorphon_classes)[3]
#' @export
rvt_geomorphon_classes <- c(
  flat = 1L, peak = 2L, ridge = 3L, shoulder = 4L, spur = 5L,
  slope = 6L, hollow = 7L, footslope = 8L, valley = 9L, pit = 10L
)

.geomorphons_tile <- function(tile, xres, yres, reach, skip, flat_threshold,
                               flat_distance, threads) {
  # skip = 0 really means "start at the adjacent cell"; the (0,0) offset is
  # dropped by .direction_offsets() anyway, so no clamp is needed here
  off <- .direction_offsets(8, reach, skip, xres, yres)
  pad <- as.integer(reach + 1L)
  geomorphons_kernel(.pad_reflect(tile, pad), pad, nrow(tile), ncol(tile),
                      off$dx, off$dy, off$dist, off$starts, off$ends,
                      flat_threshold * pi / 180, flat_distance, threads)
}

#' Geomorphons
#'
#' Classifies every cell into one of ten landform types - flat, peak, ridge,
#' shoulder, spur, slope, hollow, footslope, valley, pit - by looking at the
#' shape of the terrain around it rather than by measuring a single quantity.
#'
#' This is the odd one out in the package: the output is **categorical**, not
#' a continuous surface. Values are the class codes in
#' [rvt_geomorphon_classes], and it should be displayed with a discrete
#' palette, never a continuous ramp - the numbers are labels, and "6" is not
#' between "5" and "7" in any meaningful sense.
#'
#' What makes it useful is that the classification is *scale-explicit and
#' self-normalising*: because a cell is judged by the angles to what it can
#' see rather than by absolute heights, the same landform gets the same label
#' whether it sits in a mountain range or on a floodplain. That makes it a
#' good starting point for asking "where are all the ridges" or "where are all
#' the hollows" across ground of mixed character, and a common input to
#' further analysis rather than a picture in its own right.
#'
#' @section How it works:
#' Along each of eight directions the algorithm looks out to `reach`
#' and records the highest and lowest line-of-sight angles - the most it has to
#' look up, and the most it has to look down. Whichever is larger decides that
#' direction's verdict: terrain rises away, falls away, or is level within the
#' `flat_threshold`. That produces eight ternary votes.
#'
#' Counting how many directions rise and how many fall gives a pair of numbers
#' that indexes a lookup table of forms. Eight falling directions is a peak,
#' eight rising is a pit, half and half is a slope, and the rest of the table
#' fills in ridges, spurs, hollows, footslopes and so on. Only the counts
#' matter, not which particular directions - so the classification is
#' rotation-invariant by construction.
#'
#' It is the same ray walk as [rvt_svf()], so cost scales the same way with
#' `reach`; it just keeps both extremes per direction instead of the maximum
#' alone. Those two extremes are exactly what [rvt_openness()] and
#' [rvt_openness_negative()] average - geomorphons was introduced as an
#' extension of terrain openness, classifying the pattern of those angles
#' rather than averaging them away.
#'
#' @section Tuning:
#' * `reach` is the main control and sets **the scale of landform you get**.
#'   Small values classify micro-relief - individual banks and ditches; large
#'   values classify the hill the whole site sits on. It is in **map units**, so
#'   scale it with your resolution. There is no universally right value:
#'   choose it from the size of the landforms you care about.
#' * `flat_threshold` (degrees) decides how level counts as flat. Raise it to
#'   sweep gently undulating ground into the flat class, lower it to
#'   sub-divide subtle terrain. On noisy DEMs too small a value produces
#'   speckle.
#' * `skip` ignores the innermost ring, exactly as `noise_removal` does for
#'   [rvt_svf()], and is the right first response to a speckled result.
#' * `flat_distance` (map units, 0 to disable) relaxes the flatness threshold
#'   with range: past this distance a fixed angle would make every distant
#'   slope significant, so beyond it the threshold becomes the angle the same
#'   height difference subtends. Worth setting when `reach` is large.
#'
#' @inheritParams rvt_svf
#' @param reach how far to look, in **map units** (metres, normally;
#'   default 20)
#' @param skip innermost distance to ignore, in map units; 0 skips nothing
#'   (default 1)
#' @param flat_threshold angle below which ground counts as flat, in degrees
#'   (default 1)
#' @param flat_distance distance beyond which the flatness threshold is
#'   relaxed, in map units; 0 disables it (default 0)
#' @return `out_path`, invisibly - class codes, see [rvt_geomorphon_classes]
#' @references
#' Jasiewicz, J. and Stepinski, T. (2013) Geomorphons - a pattern recognition
#' approach to classification and mapping of landforms. *Geomorphology* 182,
#' 147-156. \doi{10.1016/j.geomorph.2012.11.005}
#' @seealso [rvt_geomorphon_classes] for the class codes; [rvt_openness()]
#'   and [rvt_openness_negative()] for the continuous measures this
#'   classifies.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' dem |> rvt_geomorphons(reach = 10)
#' @export
rvt_geomorphons <- function(dem, out_path = fs::file_temp(ext = "tif"),
                             reach = 20, skip = 1, flat_threshold = 1,
                             flat_distance = 0,
                             tile_size = NULL, threads = rvt_threads(),
                             overwrite = FALSE, progress = FALSE) {
  if (flat_threshold < 0) stop("`flat_threshold` must not be negative", call. = FALSE)
  out_path <- .out_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  reach <- .cells(reach, info$xres, "reach", "outer")
  skip <- .cells(skip, info$xres, "skip", "inner")
  if (reach < 2L)
    stop("`reach` must cover at least two cells", call. = FALSE)
  if (skip >= reach)
    stop("`skip` must be less than `reach`", call. = FALSE)

  overlap <- as.integer(reach + 1L)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

  .process_tiled(dem, list(geomorphons = out_path), overlap, tile_size,
                  function(tile, xres, yres)
                    list(geomorphons = .geomorphons_tile(tile, xres, yres, reach,
                                                          skip, flat_threshold,
                                                          flat_distance, threads)),
                  progress = progress, threads = threads,
                  # class codes: overviews by majority, never by averaging
                  categorical = TRUE)
  invisible(out_path)
}
