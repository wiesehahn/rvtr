## Slope and hillshade - both fall out of the same derivative pass, so they
## share one kernel (also used by rvt_vat(), which needs both at once).

## Passing no sun elevations makes the kernel skip the hillshade work
## entirely and return slope only.
.slope_tile <- function(tile, xres, yres, threads) {
  sa <- slope_hillshade(.pad_edge(tile, 1L), 1L, nrow(tile), ncol(tile),
                         xres, yres, numeric(0), numeric(0), threads)
  sa$slope
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
#' Slope comes from a standard finite-difference estimate over the eight
#' surrounding cells: the east-west gradient is taken from the two horizontal
#' neighbours and the north-south gradient from the two vertical ones, each
#' divided by the corresponding pixel size, and the slope is the arctangent of
#' their combined magnitude.
#'
#' Because horizontal distances come from the raster's own resolution, the
#' elevation units must match the map units (metres over metres) for the
#' angles to be right. Only the immediate neighbours are used, so this
#' responds to the finest detail - and to any noise - in the DEM. Computation
#' is at native resolution, tile by tile.
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
rvt_slope <- function(dem, out_path = tempfile(fileext = ".tif"),
                       units = c("degree", "radian", "percent"),
                       tile_size = NULL, threads = rvt_threads(),
                       overwrite = FALSE, progress = FALSE) {
  units <- match.arg(units)
  if (!overwrite && file.exists(out_path)) return(invisible(out_path))
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
rvt_hillshade <- function(dem, out_path = tempfile(fileext = ".tif"),
                           sun_azimuth = 315, sun_elevation = 35,
                           tile_size = NULL, threads = rvt_threads(),
                           overwrite = FALSE, progress = FALSE) {
  if (!overwrite && file.exists(out_path)) return(invisible(out_path))
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

#' Multi-directional hillshade
#'
#' A stack of hillshades of the same terrain lit from evenly spaced directions
#' all round the compass - one raster band per direction, each 0-1.
#'
#' A single [rvt_hillshade()] always loses the features that happen to run
#' parallel to its light. Rendering the same ground from several directions
#' means nothing stays hidden in every band, and it lets you either flick
#' between directions or combine them (an average, a maximum, or three chosen
#' bands as red/green/blue all work) to get an image with no blind direction.
#'
#' @section How it works:
#' Slope and the direction each cell faces are computed once, then reused for
#' every sun position - so `n_directions` bands cost barely more than one
#' hillshade plus the extra writing. Band `i` is lit from azimuth
#' `(360 / n_directions) * (i - 1)`, so band 1 is always due north and the
#' rest run clockwise from there, all at the same `sun_elevation`.
#'
#' The output is a single multi-band GeoTIFF. Read one band with
#' `gdalraster`, or point [rvt_plot()] at it with `band = i`.
#'
#' @section Tuning:
#' * `n_directions` sets how many bands you get. 16 is RVT's default and is
#'   generous; 4 or 8 is usually plenty to be sure nothing is hidden, and
#'   keeps the file smaller.
#' * `sun_elevation` behaves as in [rvt_hillshade()] - lower for subtle
#'   relief, higher for a flatter, more even image.
#'
#' @section Recommended settings:
#' Kokalj and Hesse (2017) use **16 directions**, with the sun elevation
#' chosen for the terrain exactly as for a single hillshade: around 35
#' degrees generally, below 10 on very flat ground, around 45 in steep or
#' complex topography.
#'
#' They pair multi-directional hillshading with a principal components
#' analysis across the directions, which is a good way to collapse the stack
#' into one or three readable images: the bands are highly correlated, so the
#' first few components carry nearly all the relief information without any
#' single direction's blind spot. That step is outside this package - take the
#' bands into your own PCA - but it is what the stack is really for.
#'
#' @references
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
#' \doi{10.3986/9789612549848}
#'
#' @inheritParams rvt_svf
#' @param n_directions number of directions (and hence bands), spread evenly
#'   around the compass starting at north (default 16)
#' @param sun_elevation height of the sun above the horizon, in degrees,
#'   shared by every direction (default 35)
#' @return `out_path`, invisibly
#' @seealso [rvt_hillshade()] for a single direction.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' rvt_multi_hillshade(dem, n_directions = 4)
#' @export
rvt_multi_hillshade <- function(dem, out_path = tempfile(fileext = ".tif"),
                                 n_directions = 16, sun_elevation = 35,
                                 tile_size = NULL, threads = rvt_threads(),
                                 overwrite = FALSE, progress = FALSE) {
  n_directions <- as.integer(n_directions)
  if (is.na(n_directions) || n_directions < 1)
    stop("`n_directions` must be at least 1", call. = FALSE)
  if (!overwrite && file.exists(out_path)) return(invisible(out_path))
  dem <- rvt_mosaic(dem)

  azimuths <- (360 / n_directions) * (seq_len(n_directions) - 1)

  info <- .dem_info(dem)
  if (is.null(tile_size)) tile_size <- .auto_tile_size(info$nx, info$ny, 1L)

  .process_tiled(dem, list(mhs = out_path), 1L, tile_size,
                  function(tile, xres, yres)
                    list(mhs = .hillshade_tile(tile, xres, yres, azimuths,
                                                sun_elevation, threads)),
                  progress = progress, threads = threads,
                  nbands = c(mhs = n_directions))
  invisible(out_path)
}
