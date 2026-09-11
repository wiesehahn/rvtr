## Time in daylight: how long each cell is lit by direct sun.

## Sun position for every time step of every supplied date, following NOAA's
## solar position algorithm (the one behind their online calculator).
##
## Accurate to well under a tenth of a degree over any plausible date range,
## which is far finer than the 10-degree azimuth bins the horizon is stored in.
## Refraction and the sun's angular width are ignored: both matter only within
## about half a degree of the horizon, where terrain shadowing dominates
## anyway.
##
## Everything is computed in UTC. For a whole day that is not an approximation
## - the set of sun positions over a full rotation is the same whatever clock
## the site keeps - and `.solar_positions()` is only ever asked for whole days.
.solar_positions <- function(dates, step_minutes, lat, lon) {
  mins <- seq(0, 1440 - step_minutes, by = step_minutes)
  m <- rep(mins, times = length(dates))
  # as.numeric(Date) counts days from 1970-01-01, whose Julian day is 2440587.5
  jd <- rep(as.numeric(dates), each = length(mins)) + 2440587.5 + m / 1440
  tc <- (jd - 2451545.0) / 36525                      # Julian centuries J2000

  d2r <- pi / 180
  l0 <- (280.46646 + tc * (36000.76983 + tc * 0.0003032)) %% 360
  ma <- (357.52911 + tc * (35999.05029 - 0.0001537 * tc)) * d2r
  ecc <- 0.016708634 - tc * (0.000042037 + 0.0000001267 * tc)

  ctr <- sin(ma) * (1.914602 - tc * (0.004817 + 0.000014 * tc)) +
         sin(2 * ma) * (0.019993 - 0.000101 * tc) +
         sin(3 * ma) * 0.000289
  omega <- (125.04 - 1934.136 * tc) * d2r
  lambda <- (l0 + ctr - 0.00569 - 0.00478 * sin(omega)) * d2r

  eps0 <- 23 + (26 + (21.448 - tc * (46.815 + tc * (0.00059 - tc * 0.001813))) / 60) / 60
  eps <- (eps0 + 0.00256 * cos(omega)) * d2r
  decl <- asin(sin(eps) * sin(lambda))

  # Equation of time (minutes): the gap between clock-uniform time and the
  # sun's own, from the Earth's tilted, elliptical orbit.
  y <- tan(eps / 2)^2
  l0r <- l0 * d2r
  eot <- 4 / d2r * (y * sin(2 * l0r) - 2 * ecc * sin(ma) +
                     4 * ecc * y * sin(ma) * cos(2 * l0r) -
                     0.5 * y * y * sin(4 * l0r) -
                     1.25 * ecc * ecc * sin(2 * ma))

  # True solar time -> hour angle, 0 at local solar noon, 15 deg per hour
  ha <- (((m + eot + 4 * lon) %% 1440) / 4 - 180) * d2r
  latr <- lat * d2r

  cos_zen <- pmin(1, pmax(-1, sin(latr) * sin(decl) +
                                cos(latr) * cos(decl) * cos(ha)))
  altitude <- 90 - acos(cos_zen) / d2r
  # Azimuth measured from south, positive westward; shift to compass bearing.
  azimuth <- (atan2(sin(ha), cos(ha) * sin(latr) - tan(decl) * cos(latr)) / d2r +
                180) %% 360

  list(azimuth = azimuth, altitude = altitude, declination = decl / d2r)
}

## Longitude and latitude of the raster's centre, for the solar geometry.
.raster_lonlat <- function(path) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  srs <- ds$getProjectionRef()
  if (!nzchar(srs))
    stop("the raster has no coordinate reference system, so its latitude ",
          "cannot be determined - pass `lat` and `lon` explicitly",
          call. = FALSE)
  bb <- ds$bbox()
  ctr <- matrix(c(mean(bb[c(1, 3)]), mean(bb[c(2, 4)])), ncol = 2)
  ll <- gdalraster::transform_xy(ctr, srs_from = srs, srs_to = "EPSG:4326")
  list(lon = ll[1, 1], lat = ll[1, 2])
}

## The dates to integrate over, and how many days the requested period spans.
## Those differ whenever `day_step` samples rather than visits every day, and
## the span is what "total over the period" has to be scaled by. Shared with
## rvt_insolation(), which takes the same period arguments.
.daylight_dates <- function(start_day, end_day, day_step, dates) {
  if (!is.null(dates)) {
    dates <- as.Date(dates)
    if (length(dates) == 0L || anyNA(dates))
      stop("`dates` must be one or more valid dates", call. = FALSE)
    return(list(dates = dates, n_period = length(dates)))
  }
  for (nm in c("start_day", "end_day", "day_step")) {
    v <- get(nm)
    if (length(v) != 1L || !is.finite(v) || v < 1 || v != floor(v))
      stop("`", nm, "` must be a single whole number of at least 1",
            call. = FALSE)
  }
  year <- as.integer(format(Sys.Date(), "%Y"))
  jan1 <- as.Date(paste0(year, "-01-01"))
  n_year <- as.integer(as.Date(paste0(year + 1L, "-01-01")) - jan1)
  start_day <- min(as.integer(start_day), n_year)
  end_day <- min(as.integer(end_day), n_year)
  # end before start wraps through New Year, so a winter season can be asked
  # for in one go
  doy <- if (end_day >= start_day) seq(start_day, end_day)
          else c(seq(start_day, n_year), seq_len(end_day))
  list(dates = jan1 + (doy[seq(1, length(doy), by = as.integer(day_step))] - 1L),
       n_period = length(doy))
}

.daylight_tile <- function(tile, ctx, xres, yres, offs, pad,
                            sun_azimuth, sun_altitude, threads) {
  padded <- .pad_reflect(tile, pad)
  do.call(daylight_kernel,
           c(list(padded, pad, nrow(tile), ncol(tile), xres, yres),
             .level_args(ctx, offs),
             list(sun_azimuth, sun_altitude, threads)))
}

#' Time in daylight
#'
#' How long each cell is lit by direct sun - the hours it actually receives,
#' after both the surrounding terrain and its own slope have taken their share.
#' A north-facing hollow in a valley may get a fraction of what the ridge above
#' it does on the same day.
#'
#' Unlike [rvt_shadow()], which freezes one instant, this integrates over time:
#' the sun is stepped across the sky and the lit steps are counted. It is an
#' analysis product rather than a picture - the input to questions about frost
#' pockets, snow lie, growing conditions, solar siting or where a site would
#' have caught the morning sun - though it also reads well as a relief image,
#' since it responds to both the shape of a slope and what stands around it.
#'
#' @section How it works:
#' Whether a cell is lit reduces to a single comparison,
#' `sun altitude > horizon(sun azimuth)`, in which the horizon depends only on
#' the compass bearing and never on the time of day. So the expensive part is
#' done once per cell and every sun position afterwards is nearly free:
#'
#' * The terrain is scanned outwards along `num_directions` compass bearings
#'   (36 by default, one every 10 degrees) out to `reach`, recording
#'   the steepest angle up to the skyline along each - the same search
#'   [rvt_svf()] uses.
#' * The cell's own slope and aspect give a second horizon: a surface tilted
#'   away from the sun goes dark well before any hill blocks it.
#' * Sun positions are computed with NOAA's solar position algorithm for every
#'   `time_step` minutes of every date in the range, from the latitude and
#'   longitude of the raster's centre. Positions below the horizon are
#'   dropped; the rest are tested against the stored horizon, interpolated
#'   between the two bearings that bracket the sun.
#'
#' Because the horizon is stored once and reused, **a finer time step or more
#' dates costs almost nothing** - the search radius is what governs the
#' runtime. Refraction and the sun's angular width are ignored, as is any
#' variation in solar geometry across the raster (one centre position is used
#' throughout), which is negligible for anything short of a continental
#' mosaic.
#'
#' @section Direct sun versus diffuse sky:
#' Three functions here measure light reaching the ground, and they answer
#' different questions:
#'
#' \describe{
#'   \item{[rvt_daylight()]}{**Direct** sunlight, **integrated over time**: how
#'     many hours the sun itself reaches a cell across a day, a season or a
#'     year. Needs a date range and a latitude, and answers in hours. A
#'     clear-sky quantity.}
#'   \item{[rvt_sky_illumination()]}{**Diffuse** skylight, at **no particular
#'     time**: how much light a cell receives from the sky dome as a whole,
#'     given both its horizon and the way it tilts. No date, no sun position,
#'     no latitude; the answer is relative, with flat open ground at 1. An
#'     overcast-sky quantity.}
#'   \item{[rvt_shadow()]}{Direct sunlight at **one instant** - a single frame
#'     of what [rvt_daylight()] integrates.}
#'   \item{[rvt_insolation()]}{The **energy** rather than the duration, in
#'     kWh/m2, with beam and diffuse together. An hour of December sun is not
#'     an hour of June sun, so this is what [rvt_daylight()]'s hour counts are
#'     really a proxy for.}
#' }
#'
#' They are complements rather than alternatives. Real ground receives both, so
#' a north-facing slope that [rvt_daylight()] shows getting no sun at all is
#' still lit diffusely, and a cell's total light is neither number by itself.
#' All three share the same horizon search, so `reach`, `num_directions` and
#' the multi-resolution levels of [rvt_reach] mean the same thing in each.
#'
#' @section Tuning:
#' * `start_day` and `end_day` choose the period, as days of the year, so the
#'   default (1 to 365) is a full year. An `end_day` earlier than `start_day`
#'   wraps through New Year, so `start_day = 305, end_day = 59` is November to
#'   February. For a single day, or a set of particular days, pass `dates`
#'   instead - `dates = as.Date("2025-06-21")` for midsummer.
#' * `day_step` samples every nth day of that period rather than all of them,
#'   and 5 is close enough to exact to treat as such: against every single day
#'   of a year, the annual mean differs by 0.006 h - under half a minute a day
#'   - for a quarter of the work. Do not raise it far. Sampling monthly looks
#'   reasonable and is not: it aliases the sun's yearly cycle and understates
#'   the annual mean by 0.3 h a day.
#' * `units` decides what the numbers mean. `"mean_hours"` (the default) is
#'   hours of direct sun on an average day of the period - the one to compare
#'   between sites or seasons. `"total_hours"` is the sum over the whole
#'   period instead, which is what you want for a seasonal energy or
#'   growing-degree style question: a full year comes to a few thousand hours.
#'   `"fraction"` divides by the daylight actually available at that latitude
#'   and season, so **1 means never shadowed** - flat open ground reads 1
#'   whatever the season.
#'
#'   `"total_hours"` scales the sampled days up to the full period, so
#'   `day_step` changes how long the calculation takes and not what the number
#'   means.
#' * `reach` is how far the terrain is allowed to cast a shadow, in map units,
#'   and is the one parameter that costs real time. It needs to cover whatever
#'   actually shades the site: a ridge 500 m away needs `reach = 500`, not the
#'   default 100. Too small and deep valleys come out too sunny.
#' * `time_step` is the step within each day, in minutes. 10 is ample; finer
#'   mostly just sharpens the quantisation of the answer.
#' * `num_directions` sets how finely the skyline is resolved. 36 matches
#'   WhiteboxTools; raising it sharpens shadow edges cast by narrow features
#'   and costs proportionally.
#'
#' @section Recommended settings:
#' For terrain analysis the annual default with `units = "mean_hours"` is the
#' usual starting point. For a specific question, use the date that poses it -
#' the winter solstice for frost hollows and snow lie, the summer solstice for
#' peak load, an equinox for a neutral average.
#'
#' Set `reach` from the terrain, not from habit: it is the distance to the
#' skyline that actually matters, and in steep country that is often several
#' hundred metres. Too short and deep valleys come out implausibly sunny.
#'
#' @inheritParams rvt_svf
#' @param start_day,end_day first and last day of the year to integrate over
#'   (1-365, defaults 1 and 365 - the whole year). `end_day` before
#'   `start_day` wraps through New Year.
#' @param day_step sample every nth day of that range (default 5)
#' @param dates specific dates to use instead, as `Date` values; overrides
#'   `start_day`/`end_day`/`day_step` when given
#' @param time_step step within each day, in minutes (default 10)
#' @param units `"mean_hours"` (default; hours of direct sun on an average
#'   day of the period), `"total_hours"` (summed over the whole period) or
#'   `"fraction"` (share of the available daylight, 1 = never shadowed)
#' @param lat,lon latitude and longitude in degrees, for the solar geometry.
#'   Taken from the centre of the raster by default, which needs the raster to
#'   carry a coordinate reference system.
#' @param num_directions compass bearings the skyline is sampled along
#'   (default 36, one every 10 degrees)
#' @param reach how far terrain can cast a shadow, in **map units**. Scanned
#'   at full resolution for the first `pyramid_px` cells (100 by default, so
#'   100 m on a 1 m DEM); past that the search reads from progressively
#'   coarser copies of the DEM, so a ridge a kilometre off is still
#'   affordable - see [rvt_reach]
#'   (default 100)
#' @return `out_path`, invisibly
#' @seealso [rvt_insolation()] for the energy rather than the duration,
#'   [rvt_shadow()] for a single instant, [rvt_sky_illumination()] for diffuse
#'   rather than direct light, and [rvt_svf()] for the horizon search they
#'   all share.
#' @examples
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#'
#' # midsummer, in hours of direct sun
#' rvt_daylight(dem, dates = as.Date("2025-06-21"), reach = 40)
#'
#' # the growing season, as a share of the daylight available
#' rvt_daylight(dem, start_day = 91, end_day = 273, units = "fraction",
#'              reach = 40)
#'
#' # total hours of sun over the same season
#' rvt_daylight(dem, start_day = 91, end_day = 273, units = "total_hours",
#'              reach = 40)
#' @export
rvt_daylight <- function(dem, out_path = tempfile(fileext = ".tif"),
                          start_day = 1, end_day = 365, day_step = 5,
                          dates = NULL, time_step = 10,
                          units = c("mean_hours", "total_hours", "fraction"),
                          lat = NULL, lon = NULL,
                          num_directions = 36, reach = 100,
                          pyramid_px = 100, pyramid_factor = 4,
                          tile_size = NULL, threads = rvt_threads(),
                          overwrite = FALSE, progress = FALSE) {
  units <- match.arg(units)
  if (!overwrite && file.exists(out_path)) return(invisible(out_path))
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

  sun_az <- sun$azimuth[up] * pi / 180
  sun_alt <- sun$altitude[up] * pi / 180

  # Counts of lit steps come back from the kernel; the scale factor turns them
  # into the requested units and is the same for every cell. "total_hours"
  # scales the sampled days up to the whole period, so `day_step` changes how
  # long it takes rather than what it means.
  per_day <- (time_step / 60) / length(dates)
  scale <- switch(units,
                   mean_hours  = per_day,
                   total_hours = per_day * n_period,
                   fraction    = 1 / sum(up))

  info <- .dem_info(dem)
  # directions are indexed by compass bearing so the kernel can map a sun
  # azimuth onto one by arithmetic
  py <- .pyramid_setup(dem, info, reach, pyramid_px, pyramid_factor, threads,
                        function(rmax, rmin, res)
                          .azimuth_offsets(num_directions, rmax, res, res,
                                            rmin = rmin))
  if (is.null(tile_size))
    tile_size <- .auto_tile_size(info$nx, info$ny, py$overlap)
  pad <- py$overlap

  .process_tiled(dem, list(daylight = out_path), py$overlap, tile_size,
                  function(tile, xres, yres, ctx = .no_aux)
                    list(daylight = .daylight_tile(tile, ctx, xres, yres,
                                                    py$offs, pad, sun_az,
                                                    sun_alt, threads) * scale),
                  progress = progress, threads = threads, aux = py$aux)
  invisible(out_path)
}
