## Slope and hillshade - both fall out of the same derivative pass, so they
## share one kernel (also used by rvt_vat(), which needs both at once).

## Passing no sun elevations makes the kernel skip the hillshade work
## entirely and return slope only.
.slope_tile <- function(tile, xres, yres, threads) {
  sa <- slope_hillshade(.pad_edge(tile, 1L), 1L, nrow(tile), ncol(tile),
                         xres, yres, numeric(0), numeric(0), threads)
  sa$slope
}

.aspect_tile <- function(tile, xres, yres, threads) {
  sa <- slope_hillshade(.pad_edge(tile, 1L), 1L, nrow(tile), ncol(tile),
                         xres, yres, numeric(0), numeric(0), threads)
  sa$aspect
}

## One hillshade per (azimuth, elevation) pair, all sharing a single pass over
## the derivatives.
.hillshade_tile <- function(tile, xres, yres, sun_azimuth, sun_elevation, threads) {
  n <- max(length(sun_azimuth), length(sun_elevation))
  sa <- slope_hillshade(.pad_edge(tile, 1L), 1L, nrow(tile), ncol(tile),
                         xres, yres, rep_len(sun_azimuth, n),
                         rep_len(sun_elevation, n), threads)
  sa$hillshade
}

#' Slope
#'
#' Steepness of the terrain at each cell, ignoring which way it faces. Flat
#' ground is 0 and a vertical face is 90 degrees (or 1.571 radians, or an
#' unbounded percentage).
#'
#' On its own this is a plain terrain derivative rather than a relief
#' visualization, but it is a useful layer to have: it picks out breaks of
#' slope very sharply, so the edges of banks, terraces, ditches and quarry
#' faces stand out as bright lines even where the features themselves are
#' low. It is also one of the ingredients [rvt_vat()] blends.
#'
#' @section How it works:
#' Slope comes from a central difference over the **four** orthogonal
#' neighbours: the east-west gradient is taken from the two horizontal ones
#' and the north-south gradient from the two vertical ones, each divided by
#' the corresponding pixel size, and the slope is the arctangent of their
#' combined magnitude.
#'
#' Several estimators are in common use and they do not give the same answer,
#' so it is worth naming the one here. This four-neighbour form responds to
#' the finest detail in the DEM, and to the most noise; Horn's eight-neighbour
#' estimate (what `gdaldem` and [gdalraster::dem_proc()] use) and the
#' Evans-Young quadratic fit average over a wider stencil and come out
#' smoother. The four-neighbour form is used here so that slope, aspect and
#' hillshade stay consistent with each other and with [rvt_vat()].
#'
#' Because horizontal distances come from the raster's own resolution, the
#' elevation units must match the map units (metres over metres) for the
#' angles to be right. Computation is at native resolution, tile by tile.
#'
#' @section Tuning:
#' There is nothing to tune but the output units. Slope has no radius or
#' smoothing parameter: it is always measured across single pixels, which
#' means a noisy DEM gives a noisy slope. If that is a problem, smooth the DEM
#' before calling this, or reach for a metric with a built-in scale such as
#' [rvt_msrm()] or [rvt_local_dominance()].
#'
#' @section Recommended settings:
#' Kokalj and Hesse (2017) recommend displaying slope in **inverted**
#' greyscale - steep dark, flat white - which keeps a surprisingly plastic,
#' relief-like impression of the terrain. Their suggested linear stretches are
#' 0-50 degrees generally, 0-15 on very flat ground, 0-60 on steep or complex
#' ground. They note Challis et al. found slope the single best technique
#' across most of the situations they tested.
#'
#' Its one serious limitation is worth stating plainly: **slope cannot tell a
#' bank from a ditch.** Rising and falling ground of the same gradient are
#' drawn identically, so a profile across a feature can easily be read as
#' convex when it is really concave. Pair it with something that carries sign
#' - [rvt_openness()] and [rvt_openness_negative()], or [rvt_msrm()], whose
#' values are signed by construction - before committing to an interpretation.
#'
#' @references
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#'
#' Pike, R. J., Evans, I. S. and Hengl, T. (2009) Geomorphometry: a brief
#' guide. In: Hengl, T. and Reuter, H. I. (eds.) *Geomorphometry: Concepts,
#' Software, Applications*. Developments in Soil Science 33. Amsterdam:
#' Elsevier, 3-30. \doi{10.1016/S0166-2481(08)00001-9}
#'
#' @inheritParams rvt_svf
#' @param units `"degree"` (default), `"radian"`, or `"percent"` (100 times
#'   the rise over the run, so 45 degrees is 100 percent)
#' @return `out_path`, invisibly
#' @seealso [rvt_hillshade()], which combines slope with the direction the
#'   terrain faces.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_slope(dem)
#' @export
rvt_slope <- function(dem, out_path = fs::file_temp(ext = "tif"),
                       units = c("degree", "radian", "percent"),
                       tile_size = NULL, threads = rvt_threads(),
                       overwrite = FALSE, progress = FALSE) {
  units <- match.arg(units)
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, 1L)

  .process_tiled(dem, list(slope = out_path), 1L, tile_size,
                  function(tile, xres, yres) {
                    slp <- .slope_tile(tile, xres, yres, threads)
                    list(slope = switch(units,
                                        radian  = slp,
                                        degree  = slp * 180 / pi,
                                        percent = tan(slp) * 100))
                  },
                  progress = progress, threads = threads)
  invisible(out_path)
}

#' Hillshade
#'
#' The familiar shaded-relief image: terrain lit by a single low sun, bright
#' where slopes face the light and dark where they turn away, from 0 to 1.
#'
#' It is the most intuitive way to look at terrain, which is why it remains
#' the default in most GIS software - but it has one well-known weakness worth
#' keeping in mind: features running parallel to the light direction almost
#' vanish, while those across it are exaggerated. If that matters, either
#' render a few different `sun_azimuth` values and compare, or use a
#' direction-free measure such as [rvt_svf()], or [rvt_asvf()] which keeps a
#' sense of direction without letting anything disappear entirely.
#'
#' @section How it works:
#' The slope and the direction each cell faces are estimated from its
#' immediate neighbours (as in [rvt_slope()]), then combined with the sun
#' position into the cosine of the angle between the surface and the incoming
#' light. Values that come out negative - surfaces pointing away from the sun
#' - are clamped to 0.
#'
#' This is local shading only: it says whether a surface is tilted towards the
#' light, not whether anything stands between it and the sun, so no long
#' shadows are cast across the terrain. Computation is at native resolution,
#' tile by tile.
#'
#' @section Tuning:
#' * `sun_azimuth` is a compass azimuth in degrees - 0 north, 90 east, 180
#'   south, 270 west. The 315 (north-west) default is a long-standing
#'   convention: lighting from the upper left is what makes relief read as
#'   raised rather than sunken to most viewers. Rotating it to run across the
#'   grain of the features you care about brings them out best.
#' * `sun_elevation` is the sun's height above the horizon in degrees. Lower
#'   values (15-25) exaggerate faint relief and are the usual choice for
#'   subtle archaeological features; higher values give a flatter, more
#'   even image.
#'
#' @section Recommended settings:
#' Kokalj and Hesse (2017) tie the sun elevation to the terrain: about
#' **35 degrees** generally, **below 10** to draw out low relief on flat or
#' gently sloping ground, and **above 45** in steep or complex topography,
#' where a low sun simply burns out one side of every hill and blacks out the
#' other. Pushed towards vertical on steep ground the image starts to resemble
#' [rvt_slope()], which they suggest as the better tool there. Azimuth stays
#' at 315 throughout; display with a linear stretch and a 2% cut-off.
#'
#' They recommend starting any investigation with a shaded relief overview,
#' because it gives the most natural-looking read of the topography and helps
#' you judge which other techniques are worth applying - but also warn about
#' three failure modes: features running parallel to the light nearly vanish
#' (so vary the azimuth, or use [rvt_multi_hillshade()]), slopes facing into
#' or away from the light saturate to white or black, and relief can read
#' *inverted*, hollows appearing as mounds. That last illusion is why lighting
#' from the north-west is conventional: most viewers read top-left lighting
#' correctly, and other azimuths invite the illusion.
#'
#' @references
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#'
#' @inheritParams rvt_svf
#' @param sun_azimuth compass azimuth of the sun, in degrees (default 315,
#'   i.e. north-west)
#' @param sun_elevation height of the sun above the horizon, in degrees
#'   (default 35)
#' @return `out_path`, invisibly
#' @seealso [rvt_asvf()] for a directional image where nothing is lost in
#'   shadow, [rvt_slope()] for the underlying steepness.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_hillshade(dem)
#' @export
rvt_hillshade <- function(dem, out_path = fs::file_temp(ext = "tif"),
                           sun_azimuth = 315, sun_elevation = 35,
                           tile_size = NULL, threads = rvt_threads(),
                           overwrite = FALSE, progress = FALSE) {
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, 1L)

  .process_tiled(dem, list(hillshade = out_path), 1L, tile_size,
                  function(tile, xres, yres)
                    list(hillshade = .hillshade_tile(tile, xres, yres,
                                                      sun_azimuth,
                                                      sun_elevation, threads)[[1]]),
                  progress = progress, threads = threads)
  invisible(out_path)
}

#' Aspect
#'
#' The compass direction each cell faces - the way water would run off it.
#' 0 is north, 90 east, 180 south, 270 west.
#'
#' Aspect on its own is rarely a useful picture: it wraps around at north, so
#' a smooth north-facing slope shows a hard seam between 359 and 0, and flat
#' ground points in an essentially arbitrary direction. It earns its place as
#' an input to other things - solar or ecological modelling, or splitting a
#' slope map into the sides of a feature.
#'
#' @section How it works:
#' Taken from the same finite-difference derivatives as [rvt_slope()], as the
#' azimuth of steepest descent, then wrapped into 0-360. Cells on perfectly
#' flat ground have no meaningful aspect; the underlying formula carries a
#' small numerical guard so they come out as a fixed direction rather than
#' undefined, which is worth knowing before you read anything into large
#' uniform patches.
#'
#' @section Tuning:
#' Nothing to tune. If you want a *readable* directional image rather than raw
#' azimuths, [rvt_hillshade()] or [rvt_asvf()] are what you actually want -
#' both fold aspect into a shaded result without the wrap-around discontinuity.
#'
#' Note this will not match `gdalraster::dem_proc(mode = "aspect")` exactly.
#' That uses Horn's eight-neighbour estimate; this uses the four-neighbour
#' central difference, so that slope and aspect here stay consistent with each
#' other and with [rvt_vat()]. On the bundled tile the two agree to about 1.4
#' degrees on slopes above 20 degrees but diverge to 33 degrees on ground
#' below 1 degree - which is not a disagreement about the terrain so much as
#' aspect being ill-defined when there is barely a gradient to point along.
#'
#' @inheritParams rvt_svf
#' @param units `"degree"` (default) or `"radian"`
#' @return `out_path`, invisibly
#' @seealso [rvt_slope()], which shares the same derivative pass.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_aspect(dem)
#' @export
rvt_aspect <- function(dem, out_path = fs::file_temp(ext = "tif"),
                        units = c("degree", "radian"),
                        tile_size = NULL, threads = rvt_threads(),
                        overwrite = FALSE, progress = FALSE) {
  units <- match.arg(units)
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  info <- .dem_info(dem)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, 1L)

  .process_tiled(dem, list(aspect = out_path), 1L, tile_size,
                  function(tile, xres, yres) {
                    a <- .aspect_tile(tile, xres, yres, threads)
                    list(aspect = if (units == "degree") a * 180 / pi else a)
                  },
                  progress = progress, threads = threads)
  invisible(out_path)
}

#' Multi-directional hillshade
#'
#' One shaded-relief image lit from several directions at once, on a 0-1
#' scale. A single [rvt_hillshade()] always loses whatever happens to run
#' parallel to its light, and buries one side of every feature in deep shadow.
#' Lighting from several sides at the same time keeps both: no direction is
#' blind, and the dark side stays readable.
#'
#' The trade is contrast: averaging several lights fills the shadows in, so
#' the image is more even but flatter than a single hillshade, and typically
#' occupies a narrower range than 0-1. A percentile stretch when plotting
#' (`rvt_plot(p, minmax_pct_cut = c(2, 98))`) gets that contrast back. Use the
#' single-direction [rvt_hillshade()] when you need a *known* light direction
#' - for a figure whose caption states it, or as the base layer of
#' [rvt_vat()].
#'
#' @section How it works:
#' A wrapper around GDAL's own `gdaldem hillshade -multidirectional`
#' (via [gdalraster::dem_proc()]), which implements the oblique-weighted
#' method of Mark (1992): the surface is lit from four sources, at azimuths
#' 225, 270, 315 and 360 degrees, all at `sun_elevation` above the horizon,
#' and the four results are averaged with weights that depend on the direction
#' each cell faces - a source lighting a slope square-on counts for more than
#' one grazing it. The single azimuth of [rvt_hillshade()] therefore has no
#' equivalent here; there is nothing to set.
#'
#' GDAL computes this in one pass over 8-neighbour (Horn) derivatives, so it
#' is fast - roughly as quick as a single hillshade - and the result is
#' rescaled here from GDAL's 1-255 bytes to a single Float32 band in 0-1, to
#' match every other product in the package. Areas of NoData in the DEM stay
#' NoData.
#'
#' @section Tuning:
#' Only `sun_elevation` and `z_factor` matter.
#'
#' * `sun_elevation` behaves as in [rvt_hillshade()]: lower picks out subtle
#'   relief more strongly, higher gives a flatter, more even image. 45 is
#'   GDAL's default and a good starting point - the weighting was designed
#'   around it - but 30-35 is worth trying when the result looks washed out.
#'   Below about 20 the image starts to look harsh even with four light
#'   sources.
#' * `z_factor` exaggerates the vertical scale. Raise it (2-5) for very flat
#'   ground; it is also the way to correct a DEM whose elevations are in
#'   different units from its coordinates.
#'
#' @section Recommended settings:
#' Kokalj and Hesse (2017) recommend multi-directional hillshading over the
#' single-direction kind for general terrain interpretation, with the sun
#' elevation chosen for the terrain much as for a single hillshade: low on
#' flat ground, around 45 degrees in steep or complex topography.
#'
#' Their own variant renders 16 separate directions and collapses them with a
#' principal components analysis. This does the same job in one band and a
#' fraction of the time; if you specifically need the PCA treatment, call
#' [rvt_hillshade()] once per azimuth and take those into your own PCA.
#'
#' @references
#' Mark, R.K. (1992) *Multidirectional, oblique-weighted, shaded-relief image
#' of the Island of Dominica*. U.S. Geological Survey Open-File Report 92-422.
#'
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#'
#' @inheritParams rvt_svf
#' @param sun_elevation height of the light sources above the horizon, in
#'   degrees, shared by all four (default 45)
#' @param z_factor vertical exaggeration applied to the elevations
#'   (default 1, no exaggeration)
#' @return `out_path`, invisibly
#' @seealso [rvt_hillshade()] for a single, named light direction.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_multi_hillshade(dem)
#' @export
rvt_multi_hillshade <- function(dem, out_path = fs::file_temp(ext = "tif"),
                                 sun_elevation = 45, z_factor = 1,
                                 tile_size = NULL, threads = rvt_threads(),
                                 overwrite = FALSE, progress = FALSE) {
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  # gdaldem tiles internally and writes Byte, with 0 reserved for NoData - no
  # valid cell ever comes out 0 - so the rescale below cannot collide with a
  # legitimately black value. -compute_edges keeps the outer ring of cells
  # instead of leaving a NoData border.
  scratch <- fs::file_temp(ext = "tif")
  on.exit(.rm_path(scratch), add = TRUE)
  gdalraster::dem_proc("hillshade", dem, scratch,
    mode_options = c("-multidirectional", "-compute_edges",
                      "-alt", format(sun_elevation, scientific = FALSE),
                      "-z", format(z_factor, scientific = FALSE)),
    quiet = TRUE)

  # Second pass only to put the result on the package's common footing:
  # Float32 in 0-1, NoData as -9999, Cloud-Optimized with overviews.
  info <- .dem_info(scratch)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, 0L)

  .process_tiled(scratch, list(mhs = out_path), 0L, tile_size,
                  function(tile, xres, yres) list(mhs = tile / 255),
                  progress = progress, threads = threads)
  invisible(out_path)
}
