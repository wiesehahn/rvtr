## Multi-scale topographic position.

.mstp_tile <- function(tile, scales, lightness, threads) {
  pad <- as.integer(max(vapply(scales, function(s) s[2], numeric(1))))
  padded <- .pad_symmetric(tile, pad)

  # RGB order is broad, meso, local - the broadest scale drives red
  bands <- lapply(scales[c("broad", "meso", "local")], function(s) {
    dev <- max_deviation_kernel(padded, pad, nrow(tile), ncol(tile),
                                 as.integer(s[1]), as.integer(s[2]),
                                 as.integer(s[3]), threads)
    v <- 1 - exp(-lightness * abs(dev))
    v[v < 0] <- 0
    v[v > 1] <- 1
    v
  })
  list(mstp = unname(bands))
}

#' Multi-scale topographic position (MSTP)
#'
#' A three-band colour composite showing, for each cell, the scale at which it
#' most stands out from its surroundings. Broad-scale prominence drives the
#' red channel, medium-scale green, and fine-scale blue - so the colour tells
#' you the *size* of the feature a cell belongs to, not just that one is
#' there.
#'
#' Reading it is largely a matter of colour: blue picks out small, sharp
#' things (field banks, ditches, individual mounds), green mid-sized landforms
#' (terraces, spurs, hollows), red the broad setting (hill masses, major
#' valleys). White means a cell stands out at every scale - a hilltop that is
#' also a local high - while dark areas are unremarkable at all three, the
#' flats and even slopes. Because each channel is scaled by how *variable* its
#' surroundings are, the result stays readable in both rugged and flat country
#' without adjustment.
#'
#' @section How it works:
#' The building block is deviation from mean elevation: for a square window,
#' `(cell - mean) / sd`. Dividing by the standard deviation is what makes the
#' three channels comparable - it measures prominence relative to local
#' roughness rather than in metres.
#'
#' For each of the three scale bands, that is evaluated at every window radius
#' from the band's minimum to its maximum, and the value with the largest
#' magnitude is kept. Each channel is then mapped to 0-1 through
#' `1 - exp(-lightness * |deviation|)`, which compresses extremes into the
#' visible range, and the three become the red, green and blue bands of a
#' single 3-band raster - so [rvt_plot()] draws it in colour by default.
#'
#' Means and standard deviations come from summed-area tables, so a
#' 2000-pixel window costs no more than a 3-pixel one; the run time is set by
#' how many radii are visited, not how large they are. Near the raster edge
#' the terrain is mirrored rather than smeared outwards - with windows this
#' large, repeating the edge value would flatten the standard deviation the
#' result divides by.
#'
#' @section Tuning:
#' * The three scales are `c(min_radius, max_radius, step)` in **map units**
#'   (metres, normally), so they mean the same ground distances whatever the
#'   resolution. The broad default reaches 2023 m, so on a small
#'   raster that channel mostly reflects the mirrored edges - crop-to-fit
#'   rather than trusting red on a tile smaller than a couple of kilometres
#'   across.
#' * `step` is purely a speed/precision trade: a coarser step visits fewer
#'   radii and runs proportionally faster, at the risk of stepping over the
#'   scale at which something stands out.
#' * `lightness` sets how quickly the colour saturates. Raise it for a
#'   brighter, higher-contrast image where moderate deviations already read
#'   strongly; lower it to reserve bright colour for the genuinely prominent.
#'
#' @inheritParams rvt_svf
#' @param local,meso,broad the three scale bands, each
#'   `c(min_radius, max_radius, step)` in **map units** (metres, normally)
#' @param lightness contrast of the colour mapping (default 1.2)
#' @return `out_path`, invisibly - a 3-band (RGB) raster
#' @seealso [rvt_msrm()], which also works across scales but returns a single
#'   band in metres rather than a colour composite.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' # scales cut down to suit a 1000 x 1000 sample tile
#' dem |> rvt_mstp(local = c(3, 21, 2), meso = c(23, 103, 18),
#'                 broad = c(123, 223, 50))
#' @export
rvt_mstp <- function(dem, out_path = fs::file_temp(ext = "tif"),
                      local = c(3, 21, 2), meso = c(23, 203, 18),
                      broad = c(223, 2023, 180), lightness = 1.2,
                      tile_size = NULL, threads = rvt_threads(),
                      overwrite = FALSE, progress = FALSE) {
  scales <- list(local = local, meso = meso, broad = broad)
  for (nm in names(scales)) {
    s <- scales[[nm]]
    if (length(s) != 3L || anyNA(s))
      stop("`", nm, "` must be c(min_radius, max_radius, step)", call. = FALSE)
    if (s[2] < s[1])
      stop("`", nm, "`: min_radius must not exceed max_radius", call. = FALSE)
  }
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  scales <- lapply(stats::setNames(names(scales), names(scales)), function(nm) {
    s <- scales[[nm]]
    c(max(1L, .cells(s[1], info$xres, paste0(nm, "[min_radius]"), "inner")),
      .cells(s[2], info$xres, paste0(nm, "[max_radius]"), "outer"),
      .cells(s[3], info$xres, paste0(nm, "[step]"), "step"))
  })
  overlap <- as.integer(max(vapply(scales, function(s) s[2], numeric(1))))
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

  .process_tiled(dem, list(mstp = out_path), overlap, tile_size,
                  function(tile, xres, yres)
                    .mstp_tile(tile, scales, lightness, threads),
                  progress = progress, threads = threads,
                  nbands = c(mstp = 3L))
  invisible(out_path)
}
