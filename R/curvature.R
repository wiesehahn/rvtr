## Curvature and Laplacian-of-Gaussian edge enhancement.

.curvature_types <- c(profile = 0L, plan = 1L, tangential = 2L,
                       mean = 3L, total = 4L, gaussian = 5L,
                       minimal = 6L, maximal = 7L, unsphericity = 8L,
                       casorati = 9L, shape_index = 10L, difference = 11L,
                       twisting = 12L, rotor = 13L)

## Two kernels, one switch: the 3x3 second differences when no radius is given
## (what ArcGIS and WhiteboxTools return), the fitted quadratic when one is.
## They are not the same at r = 1 - see src/curvature.cpp.
.curvature_tile <- function(tile, xres, yres, type, r_px, threads) {
  ty <- .curvature_types[[type]]
  if (is.null(r_px))
    return(curvature_kernel(.pad_edge(tile, 1L), 1L, nrow(tile), ncol(tile),
                             xres, yres, ty, threads))
  curvature_fit_kernel(.pad_edge(tile, r_px), r_px, nrow(tile), ncol(tile),
                        xres, yres, r_px, ty, threads)
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
#' `shape_index` is the exception: it is dimensionless and already bounded to
#' -1 to 1, so it displays as it comes.
#'
#' @section How it works:
#' A quadratic surface (Zevenbergen and Thorne 1987) is fitted to the eight
#' neighbours of each cell, and the named curvatures are combinations of its
#' second derivatives:
#'
#' \describe{
#'   \item{`profile`}{**Does water speed up or slow down here?** Positive where
#'     a slope steepens towards the brow of a bank, negative where it eases off
#'     at the foot. Measured straight down the slope. This is the one that
#'     picks out breaks of slope, and usually what people mean by "curvature"
#'     unqualified.}
#'   \item{`plan`}{**Does water gather together or spread apart here?**
#'     Negative in hollows and gullies where flow converges, positive on spurs
#'     and noses where it fans out. Measured across the slope, along the
#'     contour.}
#'   \item{`tangential`}{The same gathering-or-spreading question as `plan`,
#'     but better behaved on gentle ground, where contours wander and plan
#'     curvature turns erratic. Prefer it on anything flat.}
#'   \item{`mean`}{**Does the ground bulge outwards or dish inwards here?**
#'     Averaged over every direction at once. The general-purpose convexity
#'     measure, and the only one of these that stays meaningful on level
#'     ground.}
#'   \item{`minimal`,`maximal`}{**Which way does the ground bend most, and
#'     which way least?** A ridge bends hard across its crest and barely at all
#'     along it, so `maximal` is large and positive while `minimal` sits near
#'     zero; a valley is the mirror of that. A dome has both positive, a pit
#'     both negative.}
#'   \item{`gaussian`}{**Is this a dome or bowl, or a saddle?** Positive where
#'     the ground bends the same way in every direction, negative where it
#'     bends up one way and down the other - passes between two summits, or the
#'     junction of two valleys.}
#'   \item{`total`}{**How much is happening here at all?** How strongly the
#'     ground bends in any direction, ignoring whether it bends up or down. A
#'     magnitude, so never negative - useful for finding where something is
#'     going on when you do not care what.}
#'   \item{`unsphericity`}{**Is the bending even or lopsided?** Near zero where
#'     the ground curves the same way in every direction - a dome, a bowl, a
#'     plane - and large on ridges, valleys and saddles that bend mostly one
#'     way. Never negative.}
#'   \item{`casorati`}{**How sharp is the bending, whatever shape it makes?**
#'     Zero only where the ground is genuinely flat. Like `total` in spirit but
#'     a properly normalised curvature, so it stays comparable between gentle
#'     and steep ground.}
#'   \item{`shape_index`}{**What shape is the ground, regardless of how hard it
#'     bends?** A pit reads -1, a valley -0.5, a saddle 0, a ridge +0.5 and a
#'     peak +1. It separates form from intensity, so a faint ridge and a sharp
#'     one read alike, and being bounded it needs no stretch to display.}
#'   \item{`difference`}{**Is this place more about water speeding up, or more
#'     about water gathering?** Mean curvature averages those two and can read
#'     zero where a strong acceleration cancels a strong convergence; this is
#'     the distinction that goes missing there.}
#'   \item{`twisting`}{**Is the hillside twisting?** How fast the steepness
#'     changes as you walk *sideways* across a slope rather than up or down it.
#'     Zero on any slope of even steepness, however steep. Positive where the
#'     ground is steeper to your right as you look downhill. It is exactly the
#'     rate of change of slope angle along the contour, in radians per map
#'     unit, so the number means something directly.}
#'   \item{`rotor`}{**Does the downhill path bend as it descends?** Positive
#'     where a ball rolling down would veer right, negative left, zero where it
#'     runs straight. Picks out corkscrewing spurs and gullies that turn as
#'     they fall - a thing none of the others above can see, because they all
#'     ask about bending *across* or *along* the slope, never about the slope
#'     direction itself changing. It is the projected form of `twisting` and
#'     shares `plan`'s weakness for the same reason: on near-level ground the
#'     downhill direction is barely defined, so values explode (hundreds,
#'     against `twisting`'s fraction of one). Use `twisting` unless you
#'     specifically want the plan-view bend.}
#' }
#'
#' `unsphericity`, `casorati`, `shape_index` and `difference` come from Shary's
#' system, of which `mean`, `unsphericity` and `difference` are the three
#' independent components - every other curvature here is a combination of
#' those. `twisting` completes the basic trio of `profile`, `plan` and
#' twisting, and `rotor` is its projected form.
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
#' Give a `radius` and the scale changes: a quadratic surface is fitted by
#' least squares over the whole window (Evans 1980; Wood 1996) and the same
#' curvatures are read off its coefficients. This is the standard way to pick a
#' scale, and the only one available - the Zevenbergen and Thorne fit passes
#' exactly through nine points and has no wider form. The two are not identical
#' at the smallest radius: the fit averages the bending over every row of the
#' window, where the default takes the middle row alone. Leaving `radius` out
#' therefore keeps the values other software returns, and setting it opts into
#' the multi-scale form.
#'
#' Profile and plan curvature need a slope direction to be defined at all; on
#' genuinely flat cells they are reported as 0 rather than as a division by
#' nearly nothing. NoData inside the window contributes the centre cell's own
#' elevation, so a hole never eats a ring of cells around itself - but at a
#' large `radius` that biases curvature towards zero near big holes, so fill
#' them first with [rvt_fill()].
#'
#' @section Tuning:
#' * `type` is the only real choice, and `"profile"` is the usual starting
#'   point for finding earthworks.
#' * `radius` sets the scale, in **map units**. Omitted, curvature is measured
#'   across single pixels, which on a noisy DEM shows. A few metres steadies it
#'   and picks out banks and ditch lips; tens of metres describes the shape of
#'   the hillslope instead. Cost grows with the square of the radius, so keep
#'   it to what the features need.
#' * [rvt_log()] is the alternative when you want an edge detector rather than
#'   a curvature: its `sigma` sets a scale too, but by smoothing first.
#'
#' @inheritParams rvt_svf
#' @param type one of `"profile"` (default), `"plan"`, `"tangential"`,
#'   `"mean"`, `"minimal"`, `"maximal"`, `"gaussian"`, `"total"`,
#'   `"unsphericity"`, `"casorati"`, `"shape_index"`, `"difference"`,
#'   `"twisting"`, `"rotor"`
#' @param radius half-width of the window a quadratic surface is fitted over,
#'   in **map units**; omit it to measure across single pixels instead
#' @return `out_path`, invisibly
#' @references
#' Zevenbergen, L. W. and Thorne, C. R. (1987) Quantitative analysis of land
#' surface topography. *Earth Surface Processes and Landforms* 12, 47-56.
#'
#' Evans, I. S. (1980) An integrated system of terrain analysis and slope
#' mapping. *Zeitschrift für Geomorphologie* Supplementband 36, 274-295.
#'
#' Wood, J. (1996) *The Geomorphological Characterisation of Digital Elevation
#' Models*. PhD thesis, University of Leicester.
#'
#' Shary, P. A. (1995) Land surface in gravity points classification by a
#' complete system of curvatures. *Mathematical Geology* 27, 373-390.
#' \doi{10.1007/BF02084608}
#'
#' Koenderink, J. J. and van Doorn, A. J. (1992) Surface shape and curvature
#' scales. *Image and Vision Computing* 10, 557-564.
#' \doi{10.1016/0262-8856(92)90076-F}
#'
#' Minár, J., Evans, I. S. and Jenčo, M. (2020) A comprehensive system of
#' definitions of land surface (topographic) curvatures. *Earth-Science
#' Reviews* 211, 103414. \doi{10.1016/j.earscirev.2020.103414}
#' @seealso [rvt_log()] for an edge detector with a scale parameter.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' dem |> rvt_curvature()
#' dem |> rvt_curvature(type = "plan")
#' dem |> rvt_curvature(radius = 5)
#' @export
rvt_curvature <- function(dem, out_path = fs::file_temp(ext = "tif"),
                           type = c("profile", "plan", "tangential", "mean",
                                    "minimal", "maximal", "gaussian", "total",
                                    "unsphericity", "casorati", "shape_index",
                                    "difference", "twisting", "rotor"),
                           radius = NULL,
                           tile_size = NULL, threads = rvt_threads(),
                           overwrite = FALSE, progress = FALSE) {
  type <- match.arg(type)
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  r_px <- if (is.null(radius)) NULL
          else .cells(radius, info$xres, "radius", "outer")
  overlap <- if (is.null(r_px)) 1L else r_px
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, overlap)

  .process_tiled(dem, list(curvature = out_path), overlap, tile_size,
                  function(tile, xres, yres)
                    list(curvature = .curvature_tile(tile, xres, yres, type,
                                                      r_px, threads)),
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

## Negated: the Laplacian is positive in hollows, and this package reads convex
## as positive throughout (curvature, relief models), so LoG follows suit.
.log_tile <- function(tile, sigma, threads) {
  g <- .gauss_kernel(sigma)
  pad <- as.integer((length(g) - 1) / 2 + 1)
  -log_kernel(.pad_edge(tile, pad), pad, nrow(tile), ncol(tile), g, threads)
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
#' Values are centred on zero and **convex ground is positive**, as in
#' [rvt_curvature()] and the relief models: the lip of a bank or the top of an
#' edge reads positive, a ditch bottom or the foot of a slope negative. That is
#' the negative of the raw Laplacian, so where Kokalj and Hesse recommend an
#' inverted greyscale for the Laplacian, a plain greyscale gives the same
#' picture here. Display it with a symmetric stretch and a diverging palette.
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
#' dem |> rvt_log(sigma = 3)
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
