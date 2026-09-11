## Clear-sky solar irradiation, following the ESRA model as implemented by
## GRASS r.sun.

## Relative optical air mass at sea level (Kasten & Young 1989), from the
## refraction-corrected solar altitude. The cell's own elevation scales this
## by the air-pressure ratio, which the kernel applies.
.air_mass <- function(altitude_deg) {
  h <- altitude_deg * pi / 180
  # refraction lifts the apparent sun, most noticeably near the horizon
  dh <- 0.061359 * (0.1594 + 1.1230 * h + 0.065656 * h^2) /
        (1 + 28.9344 * h + 277.3971 * h^2)
  href <- h + dh
  1 / (sin(href) + 0.50572 * (href * 180 / pi + 6.07995)^-1.6364)
}

## Horizontal diffuse irradiance under a clear sky (ESRA): a transmission
## function of the Linke turbidity times a function of solar altitude. Both
## depend only on the sun and the atmosphere, never on the cell - which is why
## the kernel needs just their sum.
.diffuse_horizontal <- function(altitude_deg, g0, linke) {
  tn <- -0.015843 + 0.030543 * linke + 0.0003797 * linke^2
  a1p <- 0.26463 - 0.061581 * linke + 0.0031408 * linke^2
  a1 <- if (a1p * tn < 0.0022) 0.0022 / tn else a1p
  a2 <- 2.04020 + 0.018945 * linke - 0.011161 * linke^2
  a3 <- -1.3025 + 0.039231 * linke + 0.0085079 * linke^2
  s <- sin(altitude_deg * pi / 180)
  pmax(0, g0 * tn * (a1 + a2 * s + a3 * s^2))
}

## Extraterrestrial irradiance normal to the beam, corrected for the Earth's
## elliptical orbit (about +/-3.3% over the year).
.extraterrestrial <- function(dates, solar_constant = 1367) {
  doy <- as.integer(format(dates, "%j"))
  solar_constant * (1 + 0.03344 * cos(2 * pi * doy / 365.25 - 0.048869))
}

.insolation_tile <- function(tile, ctx, xres, yres, offs, pad, num_directions,
                              sun_az, sun_alt, m0, g0, linke, diffuse_sum,
                              overcast, threads) {
  do.call(insolation_kernel,
           c(list(.pad_reflect(tile, pad), pad, nrow(tile), ncol(tile),
                   xres, yres),
             .level_args(ctx, offs),
             list(sun_az, sun_alt, m0, g0,
                   .direction_azimuths(num_directions), linke, diffuse_sum,
                   overcast, threads)))
}

#' Clear-sky solar irradiation
#'
#' How much solar energy the terrain receives over a period, in kWh per square
#' metre - the quantity [rvt_daylight()]'s hour counts are a proxy for. An hour
#' of December sun is not an hour of June sun: at 47 degrees north, noon in
#' midwinter delivers about a quarter of what noon at midsummer does, and an
#' hour after sunrise is worth roughly a twelfth of noon.
#'
#' Two things cause that. Low sun crosses far more atmosphere - at 5 degrees
#' altitude the beam passes through ten atmospheres and arrives at about a
#' quarter of its strength - and it then strikes the ground at a glancing
#' angle, spreading the same beam over more surface. Both are modelled here,
#' along with the terrain's own shadowing and tilt.
#'
#' @section How it works:
#' The clear-sky model is ESRA, as implemented by GRASS GIS `r.sun`
#' (Hofierka and Suri 2002). Two components are computed and can be returned
#' separately:
#'
#' * **Beam** - direct sun. At each time step the cell is lit only if the sun
#'   clears its terrain horizon *and* the surface faces the sun. The energy
#'   arriving is the direct-normal irradiance
#'   `G0 * exp(-0.8662 * linke * m * dR(m))` - Beer-Lambert attenuation over a
#'   Rayleigh atmosphere - times the cosine of the incidence angle on the
#'   tilted surface. The air mass `m` is Kasten and Young's, corrected for
#'   refraction and for the air pressure at **the cell's own elevation**, which
#'   is worth several percent of the beam in mountains.
#' * **Diffuse** - scattered skylight. Horizontal diffuse irradiance depends
#'   only on sun altitude and turbidity, so the per-cell part is just the share
#'   of the sky dome that cell can see: exactly what [rvt_sky_illumination()]
#'   reports, with flat open ground at 1. Diffuse is *not* gated by the shadow
#'   test - ground the sun cannot reach is still lit by the sky.
#'
#' Everything else - the date range, `time_step`, `reach` and its
#' multi-resolution levels - behaves exactly as in [rvt_daylight()], and the
#' cost is the same: one horizon search per cell, with the sun positions nearly
#' free on top.
#'
#' @section On a canopy model:
#' Run on a DSM this gives the light regime of the canopy surface and of gaps,
#' with three limits worth knowing:
#'
#' * **The canopy is opaque.** Real foliage transmits and scatters, so a closed
#'   forest floor receives a few percent of open-sky radiation rather than the
#'   zero modelled here. Read the result as the geometric light regime, a lower
#'   bound.
#' * **NoData does not block.** Holes in the surface let light through, which
#'   overestimates it wherever those holes are really canopy.
#' * **Normals are noisy at crown edges**, where the 3x3 derivative sees a
#'   near-vertical facet. Gap floors are flat and unaffected.
#'
#' A gap also has to be wider than intuition suggests: the sun reaches at most
#' `90 - latitude + 23.4` degrees, so a circular gap needs a radius of about
#' `height / tan()` of that just to catch midsummer noon - at 53 degrees north,
#' 0.57 times the canopy height, and below that the direct beam is exactly zero
#' all year.
#'
#' @section This is potential, not actual, radiation:
#' There are no clouds in this model. It answers "how much does this terrain
#' allow", not "how much did this slope receive last July". Every tool of this
#' kind does the same unless fed measured radiation data, and the usual way to
#' use the output is as a relative surface - which slopes and hollows are
#' favoured - or multiplied by a measured clear-sky index for the region.
#'
#' Reflected radiation from surrounding slopes is not included either.
#'
#' @section Tuning:
#' * `linke` is the Linke turbidity factor, the atmosphere's haziness at air
#'   mass 2, and the one atmospheric parameter. Roughly: **2** very clear
#'   mountain air, **3** rural to small-town (the default, and what VOSTOK
#'   hardcodes), **4-6** urban or humid summer air. It is tabulated by region
#'   and month in the European Solar Radiation Atlas, so a defensible value can
#'   be looked up rather than guessed.
#' * `component` picks what comes back. `"global"` (default) is beam plus
#'   diffuse. `"beam"` and `"diffuse"` return one of them, and `"all"` returns
#'   a three-band raster in that order - useful because the split is where the
#'   terrain story is: beam collapses in a shaded gorge while diffuse only
#'   falls by the sky it loses.
#' * `sky_model` sets how the diffuse is spread over the sky dome. `"uniform"`
#'   (default) is the isotropic assumption normally used with a clear-sky
#'   model; `"overcast"` weights the zenith more heavily.
#' * `reach` matters more here than for a picture-making metric: a ridge that
#'   takes an hour off the winter sun changes the answer materially, so set it
#'   to the distance that actually shades the site.
#'
#' @references
#' Hofierka, J. and Suri, M. (2002) The solar radiation model for Open source
#' GIS: implementation and applications. *Proceedings of the Open Source GIS -
#' GRASS Users Conference 2002*.
#'
#' Rigollier, C., Bauer, O. and Wald, L. (2000) On the clear sky model of the
#' ESRA - European Solar Radiation Atlas - with respect to the Heliosat method.
#' *Solar Energy* 68(1), 33-48. \doi{10.1016/S0038-092X(99)00055-9}
#'
#' Kasten, F. and Young, A.T. (1989) Revised optical air mass tables and
#' approximation formula. *Applied Optics* 28(22), 4735-4738.
#'
#' @inheritParams rvt_daylight
#' @param linke Linke turbidity factor (default 3; see Tuning)
#' @param component `"global"` (default), `"beam"`, `"diffuse"`, or `"all"`
#'   for a three-band raster of beam, diffuse and global
#' @param sky_model how diffuse light is distributed over the sky dome:
#'   `"uniform"` (default) or `"overcast"`
#' @param units `"kwh"` (default), `"mj"` or `"mean_w"` - kWh/m2 and MJ/m2 are
#'   totals over the period, `"mean_w"` the mean irradiance in W/m2 across it
#' @return `out_path`, invisibly
#' @seealso [rvt_daylight()] for the duration rather than the energy, and
#'   [rvt_sky_illumination()] for the diffuse geometry on its own.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#'
#' # midsummer, kWh per square metre
#' rvt_insolation(dem, dates = as.Date("2025-06-21"), reach = 40)
#'
#' # the growing season, beam and diffuse separately
#' rvt_insolation(dem, start_day = 91, end_day = 273, component = "all",
#'                reach = 40)
#' @export
rvt_insolation <- function(dem, out_path = fs::file_temp(ext = "tif"),
                            start_day = 1, end_day = 365, day_step = 5,
                            dates = NULL, time_step = 30,
                            linke = 3,
                            component = c("global", "beam", "diffuse", "all"),
                            sky_model = c("uniform", "overcast"),
                            units = c("kwh", "mj", "mean_w"),
                            lat = NULL, lon = NULL,
                            num_directions = 36, reach = 100,
                            pyramid_px = 100, pyramid_factor = 4,
                            tile_size = NULL, threads = rvt_threads(),
                            overwrite = FALSE, progress = FALSE) {
  component <- match.arg(component)
  sky_model <- match.arg(sky_model)
  units <- match.arg(units)
  out_path <- .as_path(out_path)
  if (!overwrite && fs::file_exists(out_path)) return(invisible(out_path))
  if (!is.numeric(linke) || length(linke) != 1L || !is.finite(linke) || linke < 1)
    stop("`linke` must be a single number of at least 1", call. = FALSE)
  dem <- rvt_mosaic(dem)

  sel <- .daylight_dates(start_day, end_day, day_step, dates)
  dates <- sel$dates
  n_period <- sel$n_period

  time_step <- as.numeric(time_step)
  if (length(time_step) != 1L || is.na(time_step) || time_step <= 0 ||
      1440 %% time_step != 0)
    stop("`time_step` must be a whole number of minutes that divides 1440",
          call. = FALSE)

  if (is.null(lat) || is.null(lon)) {
    ll <- .raster_lonlat(dem)
    if (is.null(lat)) lat <- ll$lat
    if (is.null(lon)) lon <- ll$lon
  }

  sun <- .solar_positions(dates, time_step, lat, lon)
  up <- sun$altitude > 0
  if (!any(up))
    stop(sprintf("the sun never rises at latitude %.2f on any of these dates",
                  lat), call. = FALSE)

  alt <- sun$altitude[up]
  g0 <- rep(.extraterrestrial(dates), each = 1440 / time_step)[up]

  # Irradiance (W/m2) accumulated over time steps becomes energy: one step is
  # time_step/60 hours, and the sampled days stand in for the whole period.
  hours <- (time_step / 60) * (n_period / length(dates))
  scale <- switch(units,
                   kwh    = hours / 1000,
                   mj     = hours * 3600 / 1e6,
                   mean_w = 1 / sum(up))

  info <- .dem_info(dem)
  py <- .pyramid_setup(dem, info, reach, pyramid_px, pyramid_factor, threads,
                        function(rmax, rmin, res)
                          .direction_offsets(num_directions, rmax, rmin, res, res))
  if (is.null(tile_size))
    tile_size <- .auto_tile_size(info$nx, info$ny, py$overlap)
  pad <- py$overlap

  sun_az <- sun$azimuth[up] * pi / 180
  sun_alt <- alt * pi / 180
  m0 <- .air_mass(alt)
  diffuse_sum <- sum(.diffuse_horizontal(alt, g0, linke))
  overcast <- sky_model == "overcast"

  nb <- if (component == "all") 3L else 1L
  .process_tiled(dem, list(insol = out_path), py$overlap, tile_size,
                  function(tile, xres, yres, ctx = .no_aux) {
                    r <- .insolation_tile(tile, ctx, xres, yres, py$offs, pad,
                                           num_directions, sun_az, sun_alt, m0,
                                           g0, linke, diffuse_sum, overcast,
                                           threads)
                    b <- r$beam * scale
                    d <- r$diffuse * scale
                    list(insol = switch(component,
                                         global  = b + d,
                                         beam    = b,
                                         diffuse = d,
                                         all     = list(b, d, b + d)))
                  },
                  progress = progress, threads = threads, aux = py$aux,
                  nbands = c(insol = nb))
  invisible(out_path)
}
