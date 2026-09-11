## Curvature and Laplacian-of-Gaussian edge enhancement.

.curvature_types <- c(profile = 0L, plan = 1L, tangential = 2L,
                       mean = 3L, total = 4L, gaussian = 5L,
                       minimal = 6L, maximal = 7L)

.curvature_tile <- function(tile, xres, yres, type, threads) {
  curvature_kernel(.pad_edge(tile, 1L), 1L, nrow(tile), ncol(tile),
                    xres, yres, .curvature_types[[type]], threads)
}

#' Curvature
#'
#' How sharply the ground bends, and in which sense: **positive is convex**
#' (the crest of a bank, the lip of a terrace), **negative concave** (the
#' floor of a ditch, the foot of a slope), zero where the surface is locally
#' planar - whether that plane is flat or steeply tilted.
#'
#' Curvature is the natural detector for *breaks of slope*, which is often
#' what actually defines an earthwork: a bank is a convex line beside a
#' concave one. It is also the sharpest of the metrics here, in both senses -
#' it responds to the finest detail in the DEM, and it responds just as
#' eagerly to noise. Expect to need a smoothed DEM, or to prefer
#' [rvt_log()] which folds the smoothing in.
#'
#' Units are per map unit (1/m on a metric DEM), and values are small: a
#' subtle bank might reach 0.05. There is no natural range, so displaying it
#' means picking a symmetric stretch around zero - something like -0.1 to 0.1
#' - and a diverging palette so convex and concave read as opposite.
#'
#' @section How it works:
#' A quadratic surface (Zevenbergen and Thorne 1987) is fitted to the eight
#' neighbours of each cell, and the named curvatures are combinations of its
#' second derivatives:
#'
#' \describe{
#'   \item{`profile`}{curvature along the direction of steepest slope - how
#'     the gradient changes as water would run down it. This is the one that
#'     picks out breaks of slope, and usually what people mean by "curvature"
#'     unqualified.}
#'   \item{`plan`}{curvature of the contour line through the cell - whether
#'     flow converges (negative, hollows) or diverges (positive, spurs).}
#'   \item{`tangential`}{plan curvature scaled by the slope; better behaved
#'     than plan on gentle ground, where contours wander and plan curvature
#'     becomes erratic.}
#'   \item{`mean`}{the average of the two principal curvatures - a general
#'     convexity measure that stays defined on flat ground.}
#'   \item{`minimal`,`maximal`}{the two principal curvatures themselves: the
#'     gentlest and sharpest bending at that cell, over all directions. A
#'     ridge crest has a large negative minimal and a near-zero maximal; a
#'     dome has both positive.}
#'   \item{`gaussian`}{their product. Positive where the surface is dome- or
#'     bowl-like (both curvatures the same sign), negative at saddles - passes
#'     between two summits, or the junction of two valleys.}
#'   \item{`total`}{how strongly the surface bends in any direction at all,
#'     ignoring the sense - a magnitude, so never negative.}
#' }
#'
#' Except for `total`, these are the proper differential-geometry quantities,
#' carrying the `(1 + slope^2)` denominators that turn second derivatives into
#' actual curvature. Dropping those - as some implementations do - agrees only
#' where the ground is level and overstates curvature increasingly with slope,
#' by more than a factor of two by 17 degrees. `total` is the exception on
#' purpose: it follows Evans' un-normalised definition, because that is what
#' the name means in the geomorphometry literature and what other software
#' returns for it.
#'
#' Checked against closed-form values on a tilted paraboloid and against
#' WhiteboxTools: `mean`, `gaussian`, `minimal`, `maximal` and `profile` agree
#' with both to four significant figures, and `total` matches WhiteboxTools.
#'
#' Only the immediate neighbours are used, so overlap is one pixel and the
#' cost is trivial. Profile and plan curvature need a slope direction to be
#' defined at all; on genuinely flat cells they are reported as 0 rather than
#' as a division by nearly nothing.
#'
#' @section Tuning:
#' * `type` is the only real choice, and `"profile"` is the usual starting
#'   point for finding earthworks.
#' * There is no scale parameter: curvature is always measured across single
#'   pixels. On a noisy DEM that shows, and the answer is to smooth first or
#'   use [rvt_log()], whose `sigma` sets the scale explicitly.
#'
#' @inheritParams rvt_svf
#' @param type one of `"profile"` (default), `"plan"`, `"tangential"`,
#'   `"mean"`, `"minimal"`, `"maximal"`, `"gaussian"`, `"total"`
#' @return `out_path`, invisibly
#' @references
#' Zevenbergen, L. W. and Thorne, C. R. (1987) Quantitative analysis of land
#' surface topography. *Earth Surface Processes and Landforms* 12, 47-56.
#' @seealso [rvt_log()] for an edge detector with a scale parameter.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_curvature(dem)
#' rvt_curvature(dem, type = "plan")
#' @export
rvt_curvature <- function(dem, out_path = fs::file_temp(ext = "tif"),
                           type = c("profile", "plan", "tangential", "mean",
                                    "minimal", "maximal", "gaussian", "total"),
                           tile_size = NULL, threads = rvt_threads(),
                           overwrite = FALSE, progress = FALSE) {
  type <- match.arg(type)
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, 1L)

  .process_tiled(dem, list(curvature = out_path), 1L, tile_size,
                  function(tile, xres, yres)
                    list(curvature = .curvature_tile(tile, xres, yres, type, threads)),
                  progress = progress, threads = threads)
  invisible(out_path)
}

## Gaussian weights, truncated at 3 sigma (beyond that the weights are worth
## less than a thousandth of the peak and only cost time).
.gauss_kernel <- function(sigma) {
  r <- max(1L, as.integer(ceiling(3 * sigma)))
  x <- -r:r
  w <- exp(-(x^2) / (2 * sigma^2))
  w / sum(w)
}

.log_tile <- function(tile, sigma, threads) {
  g <- .gauss_kernel(sigma)
  pad <- as.integer((length(g) - 1) / 2 + 1)
  log_kernel(.pad_edge(tile, pad), pad, nrow(tile), ncol(tile), g, threads)
}

#' Laplacian of Gaussian (edge enhancement)
#'
#' Finds edges at a chosen scale: smooth the terrain by a Gaussian of width
#' `sigma`, then take its second derivative. The result is near zero on
#' smooth ground and swings sharply positive and negative on either side of a
#' break of slope, so banks, ditch lips, terrace edges and wall lines come out
#' as crisp paired lines.
#'
#' The reason to prefer this over plain [rvt_curvature()] is the `sigma`: it
#' sets which size of edge you are looking for and suppresses everything
#' finer, so a noisy lidar DTM gives a usable image instead of static.
#'
#' Kokalj and Hesse suggest using it less as a standalone product than as an
#' **overlay**: averaged into a sky-view factor or local dominance image it
#' sharpens edges and, usefully, counteracts the saturation those metrics
#' suffer on steep slopes when the stretch has been set for gentle ground.
#' Values are centred on zero, so display it with a symmetric stretch and a
#' diverging palette, or invert it - the guide notes an inverted greyscale
#' reads best.
#'
#' @section How it works:
#' The Gaussian is separable, so smoothing is a horizontal pass followed by a
#' vertical one - cost grows with `sigma`, not with `sigma` squared. The
#' kernel is truncated at three standard deviations, past which the weights
#' are below a thousandth of the peak. A 3x3 Laplacian then takes the second
#' derivative of the smoothed surface.
#'
#' NoData cells are dropped from the smoothing and the weights renormalised,
#' so a hole blurs across rather than bleeding emptiness outwards, but the
#' hole itself stays NoData. Overlap is the Gaussian radius plus one.
#'
#' @section Tuning:
#' * `sigma` is in **map units** and is the whole story: it sets the size of edge
#'   detected. Around 1-2 picks out sharp, narrow features; 4-8 finds broader
#'   terraces and scarps while ignoring texture. Kokalj and Hesse work with a
#'   Laplacian radius of about 3 cells.
#' * Larger `sigma` costs proportionally more and, more importantly, blurs
#'   away the fine detail - it is a choice about scale, not quality.
#'
#' @inheritParams rvt_svf
#' @param sigma standard deviation of the Gaussian smoothing, in **map
#'   units** (metres, normally). Unlike the window radii elsewhere it is not
#'   snapped to whole cells - it is a continuous scale - but it cannot be
#'   smaller than one cell
#'   (default 2)
#' @return `out_path`, invisibly
#' @references
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#' @seealso [rvt_curvature()] for the unsmoothed version.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_log(dem, sigma = 3)
#' @export
rvt_log <- function(dem, out_path = fs::file_temp(ext = "tif"),
                     sigma = 2, tile_size = NULL, threads = rvt_threads(),
                     overwrite = FALSE, progress = FALSE) {
  if (!is.finite(sigma) || sigma <= 0)
    stop("`sigma` must be a positive distance in map units", call. = FALSE)
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  # sigma is a continuous scale, not a window edge, so it is *not* snapped to
  # whole cells - the kernel radius it implies is rounded up instead. Below one
  # cell the Gaussian is narrower than the grid and smooths nothing.
  sigma <- sigma / info$xres
  if (sigma < 1)
    stop(sprintf("`sigma` is smaller than one cell (%s); the finest this raster supports is %s.",
                  format(info$xres), format(info$xres)), call. = FALSE)
  overlap <- as.integer((length(.gauss_kernel(sigma)) - 1) / 2 + 1)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

  .process_tiled(dem, list(log = out_path), overlap, tile_size,
                  function(tile, xres, yres) list(log = .log_tile(tile, sigma, threads)),
                  progress = progress, threads = threads)
  invisible(out_path)
}
