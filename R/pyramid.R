## Multi-resolution horizon search: near terrain at full resolution, distant
## terrain from progressively coarser copies of the same DEM.

## Level 0 is the DEM itself, searched out to `px` pixels. Each further level
## is `factor` times coarser and picks up exactly where the previous one
## stopped, again `px` pixels wide - so every level multiplies the reach by
## `factor` while costing a constant number of samples. Reaching D metres
## therefore costs O(log D) rather than O(D).
.pyramid_plan <- function(xres, reach, px = 100, factor = 4) {
  px <- as.integer(px); factor <- as.integer(factor)
  if (is.na(px) || px < 2L) stop("`pyramid_px` must be at least 2", call. = FALSE)
  if (is.na(factor) || factor < 2L)
    stop("`pyramid_factor` must be at least 2", call. = FALSE)
  if (!is.finite(reach) || reach <= 0)
    stop("`reach` must be a positive distance in map units", call. = FALSE)

  lv <- list(); inner <- 0; res <- xres
  repeat {
    outer <- min(reach, px * res)
    lv[[length(lv) + 1L]] <- list(
      fac   = as.integer(round(res / xres)),
      res   = res,
      inner = inner,
      outer = outer,
      rmin  = max(1L, .cells(inner, res, "reach", "inner")),
      rmax  = .cells(outer, res, "reach", "outer"))
    if (outer >= reach - 1e-9) break
    inner <- outer
    res <- res * factor
  }
  lv
}

## Materialise the coarse levels. `max` downsampling is not a detail: it is
## what keeps a distant skyline intact. Averaging erases any obstruction
## narrower than a coarse cell (a 3 m wall 20 m high at 300 m falls from 3.53
## to 0.65 degrees of horizon under cubic), and nearest neighbour misses such a
## wall *entirely* (0.00 degrees), since it samples one fixed position per
## block. Max is also the NoData-safe choice, picking the fill value only when
## every cell beneath is NoData.
##
## `method = "min"` is for the inverted-terrain search behind negative
## openness, where the horizon wanted is max(-z) = -min(z). Taking max(z) and
## negating would give -max(z), which is the wrong extreme.
.pyramid_build <- function(dem, plan, threads, method = "max") {
  lapply(plan, function(L) {
    if (L$fac == 1L) return(dem)
    rvt_resample(dem, L$res, method = method, threads = threads)
  })
}

## Offsets for each level, in that level's own pixels, with distances in map
## units so the angles stay comparable across levels.
.pyramid_offsets <- function(plan, num_directions, xres, yres) {
  lapply(plan, function(L)
    .direction_offsets(num_directions, L$rmax, L$rmin, L$res, L$res * yres / xres))
}

## The `aux` descriptor .process_tiled() needs: everything except level 0,
## which it already reads itself.
.pyramid_aux <- function(paths, plan) {
  if (length(plan) < 2L) return(NULL)
  lapply(seq_along(plan)[-1], function(i)
    list(path = paths[[i]], fac = plan[[i]]$fac, overlap = plan[[i]]$rmax + 1L))
}

## What .process_tiled() hands a tile function when there are no coarse levels
## at all. The kernels take the same arguments either way and simply loop over
## zero extra levels, so a short reach needs no separate code path.
.no_aux <- list(aux = list(), x0 = 0L, y0 = 0L)

## The level arguments every leveled kernel takes, in order. Keeping this in
## one place means daylight, sky illumination, shadow and the horizon kernels
## cannot drift apart in how they are called.
.level_args <- function(ctx, offs) {
  c(list(as.integer(ctx$x0), as.integer(ctx$y0),
          lapply(ctx$aux, `[[`, "m"),
          vapply(ctx$aux, `[[`, integer(1), "fac"),
          vapply(ctx$aux, `[[`, integer(1), "cx0"),
          vapply(ctx$aux, `[[`, integer(1), "cy0"),
          lapply(offs, `[[`, "dx"), lapply(offs, `[[`, "dy"),
          lapply(offs, `[[`, "dist"), lapply(offs, `[[`, "starts"),
          lapply(offs, `[[`, "ends")))
}

## Everything a leveled metric needs to set itself up: the plan, the coarse
## rasters, the per-level offsets and the tiling overlap. `offsets` builds the
## offset set for one level, given (rmax, rmin, res).
.pyramid_setup <- function(dem, info, reach, pyramid_px, pyramid_factor,
                            threads, offsets) {
  plan <- .pyramid_plan(info$xres, reach, pyramid_px, pyramid_factor)
  list(plan = plan,
       offs = lapply(plan, function(L) offsets(L$rmax, L$rmin, L$res)),
       aux = .pyramid_aux(.pyramid_build(dem, plan, threads), plan),
       overlap = as.integer(plan[[1]]$rmax + 1L))
}

#' Reaching further with a multi-resolution search
#'
#' @description
#' The horizon-search metrics - [rvt_svf()], [rvt_asvf()], [rvt_openness()],
#' [rvt_openness_negative()], [rvt_daylight()], [rvt_shadow()] and
#' [rvt_sky_illumination()] - scan outwards a fixed distance, `reach`.
#' Scanning all of it at full resolution makes the cost grow in proportion to
#' how far you look, which is why a reach long enough to catch the ridge
#' across the valley would be unaffordable on fine data.
#'
#' So they don't scan all of it at full resolution. The full-resolution scan
#' stops after `pyramid_px` cells - 100 by default - and everything beyond
#' that is read from coarser and coarser copies of the same DEM, each
#' `pyramid_factor` times coarser than the last (4 by default) and each again
#' `pyramid_px` cells wide. Cost then grows with the *logarithm* of the reach
#' instead of in proportion to it.
#'
#' # What that means for a given DEM
#'
#' That first 100 is a count of **cells**, so how far it reaches depends on the
#' resolution: 100 m on a 1 m DEM, 25 m on a 0.25 m one, 200 m on a 2 m one.
#' Ask for a `reach` within that and the search is a single full-resolution
#' scan - exactly what it has always been, and exact. Ask for more and levels
#' are added behind the scenes:
#'
#' ```
#' 0.25 m DEM, reach   20 m  ->  0.25 m data throughout
#' 0.25 m DEM, reach  100 m  ->  0.25 m out to 25 m, then 1 m to 100 m
#' 0.25 m DEM, reach 1600 m  ->  0.25 m out to 25 m, 1 m to 100 m,
#'                               4 m to 400 m, 16 m to 1600 m
#'    1 m DEM, reach  400 m  ->  1 m out to 100 m, then 4 m to 400 m
#'    2 m DEM, reach  100 m  ->  2 m data throughout
#' ```
#'
#' Nothing else about the call changes, and the coarse copies are built,
#' used and discarded for you.
#'
#' @details
#' # Why it works
#'
#' A feature's effect on the horizon is an angle, `atan(height / distance)`, so
#' detail matters in proportion to how close it is. A 1 m bump 5 m away is
#' 11 degrees of sky; the same bump 400 m away is 0.14 degrees. Sampling the
#' far field on a coarse grid therefore costs almost nothing in accuracy, while
#' saving almost all of the work.
#'
#' Level 0 is the DEM itself, searched to `pyramid_px` pixels. Each further
#' level is `pyramid_factor` times coarser and covers from where the last one
#' stopped out to `pyramid_factor` times further, again in `pyramid_px` pixels.
#' Each cell's horizon is the steepest angle found across all levels, always
#' measured from that cell's own full-resolution elevation.
#'
#' The coarse levels are built with **maximum** resampling, not averaging.
#' Averaging erases anything narrower than a coarse cell - exactly the isolated
#' crags, walls and tree lines that make up a distant skyline.
#'
#' # What it costs in accuracy
#'
#' Measured against a full-resolution search of the same reach, on terrain with
#' 300 m of relief where the 100-400 m band is worth 0.74 degrees of openness,
#' the defaults are in error by 0.068 degrees - about 9% of the far-field
#' signal they are approximating, and far less than that signal itself. On
#' gentle terrain the error is nearer 0.003 degrees.
#'
#' `pyramid_px` is the accuracy knob, and it matters much more than
#' `pyramid_factor`: doubling it roughly halves the error, because it is what
#' decides how far *full resolution* extends. Raising `pyramid_factor` is the
#' cheaper, blunter lever.
#'
#' # When to use it
#'
#' Use `reach` when distant terrain genuinely shades or encloses the site -
#' mountains, deep valleys, quarry faces - or when you want a long reach at
#' fine resolution and cannot afford it otherwise. On flat or gently rolling
#' ground the far field contributes very little (0.0004 of sky-view factor on
#' this package's sample tile), so a `reach` short enough to stay within the
#' full-resolution scan - and therefore exact - is all you need.
#'
#' `reach` is in **map units** (metres, normally), as every distance in this
#' package is, so it means the same thing whatever the resolution of the DEM.
#' To force an exact full-resolution search at any distance, raise
#' `pyramid_px` until one level covers the whole reach.
#'
#' @name rvt_reach
#' @seealso [rvt_svf()], [rvt_openness()], [rvt_daylight()], and
#'   [rvt_resample()], which builds the coarse levels.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#'
#' # look 400 m out from a 1 m DEM without paying for a 400 px search
#' dem |> rvt_svf(reach = 400)
NULL
