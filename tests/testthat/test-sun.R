# Where the sun is and what it delivers: solar geometry, time in daylight
# (rvt_daylight) and clear-sky insolation (rvt_insolation).
#
# Counting hours is not counting energy, and two results here are
# deliberately counter-intuitive: at 51 N on 21 June a north-facing slope is
# lit longer than a south-facing one, and on 21 December it inverts
# absolutely. See CLAUDE.md before reporting either as a bug.

test_that("solar positions match known geometry", {
  # noon altitude at the solstice is 90 - latitude + declination
  s <- rvtr:::.solar_positions(as.Date("2025-06-21"), 1, 51.44, 0)
  expect_equal(max(s$altitude), 90 - 51.44 + 23.44, tolerance = 0.01)
  # the equator gets ~12 h whatever the season
  for (d in c("2025-03-20", "2025-06-21", "2025-12-21")) {
    s <- rvtr:::.solar_positions(as.Date(d), 1, 0, 0)
    expect_equal(sum(s$altitude > 0) / 60, 12, tolerance = 0.02)
  }
  # at an equinox the sun rises due east
  s <- rvtr:::.solar_positions(as.Date("2025-03-20"), 1, 51.44, 0)
  up <- which(s$altitude > 0)
  expect_equal(s$azimuth[up[1]], 90, tolerance = 0.5)
})

test_that("time in daylight is exact on open flat ground", {
  flat <- make_raster(matrix(100, 60, 60))
  day <- as.Date("2025-06-21")
  lat <- 51.44; lon <- 9.83
  expected <- sum(rvtr:::.solar_positions(day, 10, lat, lon)$altitude > 0) * 10 / 60

  # tolerance is float32 rounding: results are stored as Float32
  h <- read_band(rvt_daylight(flat, dates = day, lat = lat, lon = lon,
                               reach = 20))
  expect_equal(max(h), expected, tolerance = 1e-6)
  expect_equal(min(h), expected, tolerance = 1e-6)

  f <- read_band(rvt_daylight(flat, dates = day, lat = lat, lon = lon,
                               reach = 20, units = "fraction"))
  expect_equal(max(f), 1, tolerance = 1e-6)
})

test_that("a shaft gets no sun while the plateau around it gets all of it", {
  m <- matrix(500, 60, 60)
  m[28:32, 28:32] <- 0
  p <- make_raster(m)
  day <- as.Date("2025-06-21")
  h <- read_band(rvt_daylight(p, dates = day, lat = 51.44, lon = 9.83,
                               reach = 20))
  expect_equal(h[30, 30], 0)
  expect_gt(h[5, 5], 16)
})

test_that("slope aspect and season interact the way the sun does", {
  # A 30 degree slope at 51 N. In midsummer the sun rises in the north-east
  # and is 62 degrees high at noon, so the north-facing slope is lit all day
  # while the south-facing one is self-shadowed at dawn and dusk. Come the
  # winter solstice the sun peaks at 15 degrees and the ordering inverts
  # completely. Counting hours is not the same as counting energy.
  n <- 80
  fall <- tan(30 * pi / 180)
  south <- outer(seq(n, 1) * fall, rep(1, n))   # drops southward
  north <- outer(seq(1, n) * fall, rep(1, n))   # drops northward
  hours <- function(m, d) {
    read_band(rvt_daylight(make_raster(m), dates = as.Date(d), lat = 51.44,
                            lon = 9.83, reach = 20))[40, 40]
  }
  expect_gt(hours(north, "2025-06-21"), hours(south, "2025-06-21"))
  expect_equal(hours(north, "2025-12-21"), 0)
  expect_gt(hours(south, "2025-12-21"), 7)
})

test_that("the day range selects the days it says it does", {
  yr <- as.integer(format(Sys.Date(), "%Y"))
  jan1 <- as.Date(paste0(yr, "-01-01"))
  # a stand-in for the date building inside rvt_daylight()
  span <- function(...) {
    f <- formals(rvt_daylight)
    a <- list(...)
    n_year <- as.integer(as.Date(paste0(yr + 1L, "-01-01")) - jan1)
    s <- if (is.null(a$start_day)) f$start_day else a$start_day
    e <- if (is.null(a$end_day)) f$end_day else a$end_day
    st <- if (is.null(a$day_step)) f$day_step else a$day_step
    doy <- if (e >= s) seq(s, e) else c(seq(s, n_year), seq_len(e))
    jan1 + (doy[seq(1, length(doy), by = st)] - 1L)
  }
  expect_equal(length(span()), 73L)                  # whole year, every 5th day
  expect_equal(min(span()), jan1)
  expect_equal(length(span(start_day = 1, end_day = 10, day_step = 1)), 10L)
  # a range ending before it starts wraps through New Year
  w <- span(start_day = 360, end_day = 5, day_step = 1)
  expect_equal(length(w), 11L)
  expect_true(any(format(w, "%m") == "12") && any(format(w, "%m") == "01"))

  expect_error(rvt_daylight(dem, start_day = 0), "at least 1")
  expect_error(rvt_daylight(dem, day_step = 2.5), "whole number")
})

test_that("a narrower day range gives a different answer than the full year", {
  # cheap settings: a handful of days, coarse time step, short reach
  june <- read_band(rvt_daylight(dem, start_day = 152, end_day = 182,
                                  day_step = 15, time_step = 60, reach = 10))
  dec <- read_band(rvt_daylight(dem, start_day = 335, end_day = 365,
                                 day_step = 15, time_step = 60, reach = 10))
  expect_gt(mean(june, na.rm = TRUE), mean(dec, na.rm = TRUE))
  # midsummer must beat midwinter everywhere at this latitude
  expect_true(all(june >= dec, na.rm = TRUE))
})

test_that("total hours is the sum over the period, mean hours the daily average", {
  flat <- make_raster(matrix(100, 60, 60))
  lat <- 51.44; lon <- 9.83
  # a short, fully sampled span so the total is unambiguous
  args <- list(flat, start_day = 100, end_day = 109, day_step = 1,
               time_step = 30, lat = lat, lon = lon, reach = 20)

  mean_h <- read_band(do.call(rvt_daylight, c(args, list(units = "mean_hours"))))
  total_h <- read_band(do.call(rvt_daylight, c(args, list(units = "total_hours"))))
  expect_equal(max(total_h), max(mean_h) * 10, tolerance = 1e-5)

  # on open flat ground the total must equal the astronomical day lengths added up
  dates <- as.Date(paste0(format(Sys.Date(), "%Y"), "-01-01")) + (100:109 - 1)
  expected <- sum(rvtr:::.solar_positions(dates, 30, lat, lon)$altitude > 0) * 30 / 60
  expect_equal(max(total_h), expected, tolerance = 1e-5)

  # sampling fewer days must not change what "total" means, only its precision
  sparse <- read_band(do.call(rvt_daylight,
    c(args[names(args) != "day_step"], list(day_step = 5, units = "total_hours"))))
  expect_equal(max(sparse), max(total_h), tolerance = 0.02)
})

test_that("the clear-sky solar geometry matches published values", {
  # Kasten & Young air mass: 1 at the zenith, ~2 at 30 degrees, ~10 at 5
  expect_equal(.air_mass(90), 1, tolerance = 1e-3)
  expect_equal(.air_mass(30), 2, tolerance = 0.02)
  expect_equal(.air_mass(5), 10, tolerance = 0.1)
  # and it must grow monotonically as the sun drops
  a <- .air_mass(c(90, 60, 30, 10, 5))
  expect_true(all(diff(a) > 0))

  # Earth-Sun distance: perihelion in early January, about +/-3.3%
  e <- .extraterrestrial(as.Date(c("2025-01-03", "2025-07-04")))
  expect_gt(e[1], e[2])
  expect_equal(e[1] / e[2], 1.069, tolerance = 0.01)

  # diffuse rises with sun altitude and with haze
  expect_gt(.diffuse_horizontal(60, 1367, 3), .diffuse_horizontal(20, 1367, 3))
  expect_gt(.diffuse_horizontal(40, 1367, 6), .diffuse_horizontal(40, 1367, 2))
})

test_that("insolation on open flat ground is physically plausible", {
  flat <- make_raster(matrix(700, 40, 40))
  args <- list(flat, lat = 47.5, lon = 11.1, reach = 20, num_directions = 16)

  day <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), time_step = 10))))
  # a cloudless June day at 47 N delivers roughly 7-9 kWh/m2 on the horizontal
  expect_gt(mean(day), 7)
  expect_lt(mean(day), 10)

  # beam plus diffuse must equal global, cell for cell
  b <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), component = "beam"))))
  d <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), component = "diffuse"))))
  g <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), component = "global"))))
  expect_equal(b + d, g, tolerance = 1e-5)
  expect_gt(mean(b), mean(d))              # clear sky: beam dominates

  # "all" returns the same three as bands
  p <- do.call(rvt_insolation, c(args, list(dates = as.Date("2025-06-21"),
                                             component = "all")))
  ds <- methods::new(gdalraster::GDALRaster, p, read_only = TRUE)
  expect_equal(ds$getRasterCount(), 3L)
  ds$close()

  # haze cuts the beam and lifts the diffuse
  hazy_b <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), component = "beam", linke = 6))))
  hazy_d <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), component = "diffuse", linke = 6))))
  expect_lt(mean(hazy_b), mean(b))
  expect_gt(mean(hazy_d), mean(d))

  # units are just a rescaling of the same field
  mj <- read_band(do.call(rvt_insolation,
    c(args, list(dates = as.Date("2025-06-21"), units = "mj"))))
  expect_equal(mj, g * 3.6, tolerance = 1e-4)
})

test_that("insolation responds to aspect, and more sharply in winter", {
  n <- 60
  fall <- tan(30 * pi / 180)
  south <- outer(seq(n, 1) * fall, rep(1, n)) + 700
  north <- outer(seq(1, n) * fall, rep(1, n)) + 700
  mid <- function(m) mean(m[25:35, 25:35])
  val <- function(z, d) mid(read_band(rvt_insolation(make_raster(z),
    dates = as.Date(d), lat = 47.5, lon = 11.1, reach = 20,
    num_directions = 16)))

  s6 <- val(south, "2025-06-21"); n6 <- val(north, "2025-06-21")
  s12 <- val(south, "2025-12-21"); n12 <- val(north, "2025-12-21")

  expect_gt(s6, n6)                        # south wins in both seasons
  expect_gt(s12, n12)
  expect_gt(s12 / n12, s6 / n6)            # but far more decisively in winter

  # in midwinter the north slope gets no direct sun at all, yet is not black:
  # the sky still lights it, and diffuse is not gated by the shadow test
  nb <- mid(read_band(rvt_insolation(make_raster(north),
    dates = as.Date("2025-12-21"), component = "beam", lat = 47.5, lon = 11.1,
    reach = 20, num_directions = 16)))
  expect_lt(nb, 0.05)
  expect_gt(n12, 0.2)
})

