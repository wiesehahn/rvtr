# Time in daylight

How long each cell is lit by direct sun - the hours it actually
receives, after both the surrounding terrain and its own slope have
taken their share. A north-facing hollow in a valley may get a fraction
of what the ridge above it does on the same day.

## Usage

``` r
rvt_daylight(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  start_day = 1,
  end_day = 365,
  day_step = 5,
  dates = NULL,
  time_step = 10,
  units = c("mean_hours", "total_hours", "fraction"),
  lat = NULL,
  lon = NULL,
  num_directions = 36,
  reach = 100,
  pyramid_px = 100,
  pyramid_factor = 4,
  pyramid_method = "max",
  tile_size = NULL,
  threads = rvt_threads(),
  overwrite = FALSE,
  progress = FALSE
)
```

## Arguments

- dem:

  path to a DEM, or a vector of paths forming a mosaic (see
  [`rvt_mosaic()`](https://wiesehahn.github.io/rvtr/reference/rvt_mosaic.md));
  mosaics are read across file boundaries, so tiles do not produce
  seams. A terra `SpatRaster` works too: one read from a file is used as
  it lies, while one terra computed in memory is written to a temporary
  GeoTIFF first.

- out_path:

  output GeoTIFF path (default: a temp file). A folder that does not
  exist yet is created.

- start_day, end_day:

  first and last day of the year to integrate over (1-365, defaults 1
  and 365 - the whole year). `end_day` before `start_day` wraps through
  New Year.

- day_step:

  sample every nth day of that range (default 5)

- dates:

  specific dates to use instead, as `Date` values; overrides
  `start_day`/`end_day`/`day_step` when given

- time_step:

  step within each day, in minutes (default 10)

- units:

  `"mean_hours"` (default; hours of direct sun on an average day of the
  period), `"total_hours"` (summed over the whole period) or
  `"fraction"` (share of the available daylight, 1 = never shadowed)

- lat, lon:

  latitude and longitude in degrees, for the solar geometry. Taken from
  the centre of the raster by default, which needs the raster to carry a
  coordinate reference system.

- num_directions:

  compass bearings the skyline is sampled along (default 36, one every
  10 degrees)

- reach:

  how far terrain can cast a shadow, in **map units**. Scanned at full
  resolution for the first `pyramid_px` cells (100 by default, so 100 m
  on a 1 m DEM); past that the search reads from progressively coarser
  copies of the DEM, so a ridge a kilometre off is still affordable -
  see
  [rvt_reach](https://wiesehahn.github.io/rvtr/reference/rvt_reach.md)
  (default 100)

- pyramid_px, pyramid_factor:

  shape of that multi-resolution search: pixels scanned per level
  (default 100, the accuracy knob) and how much coarser each level is
  than the last (default 4). Only used with `reach`.

- pyramid_method:

  how the coarse levels are built. `"max"` (default) is right for
  bare-earth terrain, where a narrow wall or crag must survive
  coarsening. On a **rough surface - a canopy or vegetation model - use
  `"q3"`**: max there takes the tallest crown in each block and treats
  it as solid, over-stating obstruction, and measures about twice the
  error of `"q3"` against a full-resolution search. See
  [rvt_reach](https://wiesehahn.github.io/rvtr/reference/rvt_reach.md).

- tile_size:

  tile edge in pixels; NULL picks a size targeting roughly 256 MB per
  working matrix

- threads:

  C++ threads (see
  [`rvt_threads()`](https://wiesehahn.github.io/rvtr/reference/rvt_threads.md))

- overwrite:

  recompute even if `out_path` exists

- progress:

  report per-tile progress

## Value

`out_path`, invisibly

## Details

Unlike
[`rvt_shadow()`](https://wiesehahn.github.io/rvtr/reference/rvt_shadow.md),
which freezes one instant, this integrates over time: the sun is stepped
across the sky and the lit steps are counted. It is an analysis product
rather than a picture - the input to questions about frost pockets, snow
lie, growing conditions, solar siting or where a site would have caught
the morning sun - though it also reads well as a relief image, since it
responds to both the shape of a slope and what stands around it.

## How it works

Whether a cell is lit reduces to a single comparison,
`sun altitude > horizon(sun azimuth)`, in which the horizon depends only
on the compass bearing and never on the time of day. So the expensive
part is done once per cell and every sun position afterwards is nearly
free:

- The terrain is scanned outwards along `num_directions` compass
  bearings (36 by default, one every 10 degrees) out to `reach`,
  recording the steepest angle up to the skyline along each - the same
  search
  [`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md)
  uses.

- The cell's own slope and aspect give a second horizon: a surface
  tilted away from the sun goes dark well before any hill blocks it.

- Sun positions are computed with NOAA's solar position algorithm for
  every `time_step` minutes of every date in the range, from the
  latitude and longitude of the raster's centre. Positions below the
  horizon are dropped; the rest are tested against the stored horizon,
  interpolated between the two bearings that bracket the sun.

Because the horizon is stored once and reused, **a finer time step or
more dates costs almost nothing** - the search radius is what governs
the runtime. Refraction and the sun's angular width are ignored, as is
any variation in solar geometry across the raster (one centre position
is used throughout), which is negligible for anything short of a
continental mosaic.

## Direct sun versus diffuse sky

Three functions here measure light reaching the ground, and they answer
different questions:

- `rvt_daylight()`:

  **Direct** sunlight, **integrated over time**: how many hours the sun
  itself reaches a cell across a day, a season or a year. Needs a date
  range and a latitude, and answers in hours. A clear-sky quantity.

- [`rvt_sky_illumination()`](https://wiesehahn.github.io/rvtr/reference/rvt_sky_illumination.md):

  **Diffuse** skylight, at **no particular time**: how much light a cell
  receives from the sky dome as a whole, given both its horizon and the
  way it tilts. No date, no sun position, no latitude; the answer is
  relative, with flat open ground at 1. An overcast-sky quantity.

- [`rvt_shadow()`](https://wiesehahn.github.io/rvtr/reference/rvt_shadow.md):

  Direct sunlight at **one instant** - a single frame of what
  `rvt_daylight()` integrates.

- [`rvt_insolation()`](https://wiesehahn.github.io/rvtr/reference/rvt_insolation.md):

  The **energy** rather than the duration, in kWh/m2, with beam and
  diffuse together. An hour of December sun is not an hour of June sun,
  so this is what `rvt_daylight()`'s hour counts are really a proxy for.

They are complements rather than alternatives. Real ground receives
both, so a north-facing slope that `rvt_daylight()` shows getting no sun
at all is still lit diffusely, and a cell's total light is neither
number by itself. All three share the same horizon search, so `reach`,
`num_directions` and the multi-resolution levels of
[rvt_reach](https://wiesehahn.github.io/rvtr/reference/rvt_reach.md)
mean the same thing in each.

## Tuning

- `start_day` and `end_day` choose the period, as days of the year, so
  the default (1 to 365) is a full year. An `end_day` earlier than
  `start_day` wraps through New Year, so `start_day = 305, end_day = 59`
  is November to February. For a single day, or a set of particular
  days, pass `dates` instead - `dates = as.Date("2025-06-21")` for
  midsummer.

- `day_step` samples every nth day of that period rather than all of
  them, and 5 is close enough to exact to treat as such: against every
  single day of a year, the annual mean differs by 0.006 h - under half
  a minute a day

  - for a quarter of the work. Do not raise it far. Sampling monthly
    looks reasonable and is not: it aliases the sun's yearly cycle and
    understates the annual mean by 0.3 h a day.

- `units` decides what the numbers mean. `"mean_hours"` (the default) is
  hours of direct sun on an average day of the period - the one to
  compare between sites or seasons. `"total_hours"` is the sum over the
  whole period instead, which is what you want for a seasonal energy or
  growing-degree style question: a full year comes to a few thousand
  hours. `"fraction"` divides by the daylight actually available at that
  latitude and season, so **1 means never shadowed** - flat open ground
  reads 1 whatever the season.

  `"total_hours"` scales the sampled days up to the full period, so
  `day_step` changes how long the calculation takes and not what the
  number means.

- `reach` is how far the terrain is allowed to cast a shadow, in map
  units, and is the one parameter that costs real time. It needs to
  cover whatever actually shades the site: a ridge 500 m away needs
  `reach = 500`, not the default 100. Too small and deep valleys come
  out too sunny.

- `time_step` is the step within each day, in minutes. 10 is ample;
  finer mostly just sharpens the quantisation of the answer.

- `num_directions` sets how finely the skyline is resolved. 36 matches
  WhiteboxTools; raising it sharpens shadow edges cast by narrow
  features and costs proportionally.

## Recommended settings

For terrain analysis the annual default with `units = "mean_hours"` is
the usual starting point. For a specific question, use the date that
poses it - the winter solstice for frost hollows and snow lie, the
summer solstice for peak load, an equinox for a neutral average.

Set `reach` from the terrain, not from habit: it is the distance to the
skyline that actually matters, and in steep country that is often
several hundred metres. Too short and deep valleys come out implausibly
sunny.

## See also

[`rvt_insolation()`](https://wiesehahn.github.io/rvtr/reference/rvt_insolation.md)
for the energy rather than the duration,
[`rvt_shadow()`](https://wiesehahn.github.io/rvtr/reference/rvt_shadow.md)
for a single instant,
[`rvt_sky_illumination()`](https://wiesehahn.github.io/rvtr/reference/rvt_sky_illumination.md)
for diffuse rather than direct light, and
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) for
the horizon search they all share.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")

# midsummer, in hours of direct sun
dem |> rvt_daylight(dates = as.Date("2025-06-21"), reach = 40)

# the growing season, as a share of the daylight available
dem |> rvt_daylight(start_day = 91, end_day = 273, units = "fraction",
                    reach = 40)

# total hours of sun over the same season
dem |> rvt_daylight(start_day = 91, end_day = 273, units = "total_hours",
                    reach = 40)
```
