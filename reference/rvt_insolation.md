# Clear-sky solar irradiation

How much solar energy the terrain receives over a period, in kWh per
square metre - the quantity
[`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md)'s
hour counts are a proxy for. An hour of December sun is not an hour of
June sun: at 47 degrees north, noon in midwinter delivers about a
quarter of what noon at midsummer does, and an hour after sunrise is
worth roughly a twelfth of noon.

## Usage

``` r
rvt_insolation(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  start_day = 1,
  end_day = 365,
  day_step = 5,
  dates = NULL,
  time_step = 30,
  linke = 3,
  component = c("global", "beam", "diffuse", "all"),
  sky_model = c("uniform", "overcast"),
  units = c("kwh", "mj", "mean_w"),
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

  output GeoTIFF path (default: a temp file)

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

- linke:

  Linke turbidity factor (default 3; see Tuning)

- component:

  `"global"` (default), `"beam"`, `"diffuse"`, or `"all"` for a
  three-band raster of beam, diffuse and global

- sky_model:

  how diffuse light is distributed over the sky dome: `"uniform"`
  (default) or `"overcast"`

- units:

  `"kwh"` (default), `"mj"` or `"mean_w"` - kWh/m2 and MJ/m2 are totals
  over the period, `"mean_w"` the mean irradiance in W/m2 across it

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

Two things cause that. Low sun crosses far more atmosphere - at 5
degrees altitude the beam passes through ten atmospheres and arrives at
about a quarter of its strength - and it then strikes the ground at a
glancing angle, spreading the same beam over more surface. Both are
modelled here, along with the terrain's own shadowing and tilt.

## How it works

The clear-sky model is ESRA, as implemented by GRASS GIS `r.sun`
(Hofierka and Suri 2002). Two components are computed and can be
returned separately:

- **Beam** - direct sun. At each time step the cell is lit only if the
  sun clears its terrain horizon *and* the surface faces the sun. The
  energy arriving is the direct-normal irradiance
  `G0 * exp(-0.8662 * linke * m * dR(m))` - Beer-Lambert attenuation
  over a Rayleigh atmosphere - times the cosine of the incidence angle
  on the tilted surface. The air mass `m` is Kasten and Young's,
  corrected for refraction and for the air pressure at **the cell's own
  elevation**, which is worth several percent of the beam in mountains.

- **Diffuse** - scattered skylight. Horizontal diffuse irradiance
  depends only on sun altitude and turbidity, so the per-cell part is
  just the share of the sky dome that cell can see: exactly what
  [`rvt_sky_illumination()`](https://wiesehahn.github.io/rvtr/reference/rvt_sky_illumination.md)
  reports, with flat open ground at 1. Diffuse is *not* gated by the
  shadow test - ground the sun cannot reach is still lit by the sky.

Everything else - the date range, `time_step`, `reach` and its
multi-resolution levels - behaves exactly as in
[`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md),
and the cost is the same: one horizon search per cell, with the sun
positions nearly free on top.

## On a canopy model

Run on a DSM this gives the light regime of the canopy surface and of
gaps, with three limits worth knowing:

- **The canopy is opaque.** Real foliage transmits and scatters, so a
  closed forest floor receives a few percent of open-sky radiation
  rather than the zero modelled here. Read the result as the geometric
  light regime, a lower bound.

- **NoData does not block.** Holes in the surface let light through,
  which overestimates it wherever those holes are really canopy.

- **Normals are noisy at crown edges**, where the 3x3 derivative sees a
  near-vertical facet. Gap floors are flat and unaffected.

- **Set `pyramid_method = "q3"`** if `reach` is long enough to build
  coarse levels. The default `"max"` takes the tallest crown in each
  coarse block and treats it as solid, which over-states shading on a
  surface that is all peaks; see
  [rvt_reach](https://wiesehahn.github.io/rvtr/reference/rvt_reach.md).

A gap also has to be wider than intuition suggests: the sun reaches at
most `90 - latitude + 23.4` degrees, so a circular gap needs a radius of
about `height / tan()` of that just to catch midsummer noon - at 53
degrees north, 0.57 times the canopy height, and below that the direct
beam is exactly zero all year.

## This is potential, not actual, radiation

There are no clouds in this model. It answers "how much does this
terrain allow", not "how much did this slope receive last July". Every
tool of this kind does the same unless fed measured radiation data, and
the usual way to use the output is as a relative surface - which slopes
and hollows are favoured - or multiplied by a measured clear-sky index
for the region.

Reflected radiation from surrounding slopes is not included either.

## Tuning

- `linke` is the Linke turbidity factor, the atmosphere's haziness at
  air mass 2, and the one atmospheric parameter. Roughly: **2** very
  clear mountain air, **3** rural to small-town (the default, and what
  VOSTOK hardcodes), **4-6** urban or humid summer air. It is tabulated
  by region and month in the European Solar Radiation Atlas, so a
  defensible value can be looked up rather than guessed.

- `component` picks what comes back. `"global"` (default) is beam plus
  diffuse. `"beam"` and `"diffuse"` return one of them, and `"all"`
  returns a three-band raster in that order - useful because the split
  is where the terrain story is: beam collapses in a shaded gorge while
  diffuse only falls by the sky it loses.

- `sky_model` sets how the diffuse is spread over the sky dome.
  `"uniform"` (default) is the isotropic assumption normally used with a
  clear-sky model; `"overcast"` weights the zenith more heavily.

- `reach` matters more here than for a picture-making metric: a ridge
  that takes an hour off the winter sun changes the answer materially,
  so set it to the distance that actually shades the site.

## References

Hofierka, J. and Suri, M. (2002) The solar radiation model for Open
source GIS: implementation and applications. *Proceedings of the Open
Source GIS - GRASS Users Conference 2002*.

Rigollier, C., Bauer, O. and Wald, L. (2000) On the clear sky model of
the ESRA - European Solar Radiation Atlas - with respect to the Heliosat
method. *Solar Energy* 68(1), 33-48.
[doi:10.1016/S0038-092X(99)00055-9](https://doi.org/10.1016/S0038-092X%2899%2900055-9)

Kasten, F. and Young, A.T. (1989) Revised optical air mass tables and
approximation formula. *Applied Optics* 28(22), 4735-4738.

## See also

[`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md)
for the duration rather than the energy, and
[`rvt_sky_illumination()`](https://wiesehahn.github.io/rvtr/reference/rvt_sky_illumination.md)
for the diffuse geometry on its own.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")

# midsummer, kWh per square metre
dem |> rvt_insolation(dates = as.Date("2025-06-21"), reach = 40)

# the growing season, beam and diffuse separately
dem |> rvt_insolation(start_day = 91, end_day = 273, component = "all",
                      reach = 40)
```
