# Cast shadow

Which ground is lit and which lies in shadow for a given sun position: 1
where the sun reaches the cell, 0 where terrain between the cell and the
sun blocks it.

## Usage

``` r
rvt_shadow(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  sun_azimuth = 315,
  sun_elevation = 35,
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

- sun_azimuth:

  compass azimuth of the sun, in degrees (default 315)

- sun_elevation:

  height of the sun above the horizon, in degrees (default 35)

- reach:

  how far along the ray to look for blocking terrain, in **map units**
  (metres, normally; default 100). Scanned at full resolution for the
  first `pyramid_px` cells (100 by default, so 100 m on a 1 m DEM); past
  that the search reads from progressively coarser copies of the DEM, so
  a shadow cast from kilometres away stays affordable - see
  [rvt_reach](https://wiesehahn.github.io/rvtr/reference/rvt_reach.md)

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

`out_path`, invisibly - 1 lit, 0 shadowed

## Details

This is genuine cast shadow, not the shading in
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md).
Hillshade only asks whether a surface is tilted towards the light, so it
never puts anything in the shadow of a hill in front of it; this traces
the terrain towards the sun and finds out. Low sun angles throw long
shadows that pick out very slight relief, which is why raking light is a
standard trick for spotting earthworks.

## How it works

A single ray is traced from each cell towards `sun_azimuth`, out to
`reach`, keeping the steepest angle up to the terrain along it - the
horizon in the sun's direction. If that horizon stands higher than
`sun_elevation`, the sun is blocked and the cell is in shadow.

Only one direction is searched, so this is far cheaper than the
all-round metrics: roughly one thirty-second of
[`rvt_sky_illumination()`](https://wiesehahn.github.io/rvtr/reference/rvt_sky_illumination.md)
at the same radius. The output is 0 or 1, with NoData preserved.

## Tuning

- `reach` sets how far a shadow can reach. A shadow cast by terrain
  further away than this is missed, so it needs to cover the tallest
  feature divided by the tangent of the sun elevation - a 20 m rise at
  10 degrees casts a shadow about 113 m long. Low sun therefore needs a
  large radius to be honest.

- `sun_elevation` is the main control on the picture: low values (5-15)
  throw long shadows that reveal very subtle relief, high values shadow
  only steep ground.

- `sun_azimuth` is a compass azimuth, as in
  [`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md).

## Direct sun versus diffuse sky

Three functions here measure light reaching the ground, and they answer
different questions:

- [`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md):

  **Direct** sunlight, **integrated over time**: how many hours the sun
  itself reaches a cell across a day, a season or a year. Needs a date
  range and a latitude, and answers in hours. A clear-sky quantity.

- [`rvt_sky_illumination()`](https://wiesehahn.github.io/rvtr/reference/rvt_sky_illumination.md):

  **Diffuse** skylight, at **no particular time**: how much light a cell
  receives from the sky dome as a whole, given both its horizon and the
  way it tilts. No date, no sun position, no latitude; the answer is
  relative, with flat open ground at 1. An overcast-sky quantity.

- `rvt_shadow()`:

  Direct sunlight at **one instant** - a single frame of what
  [`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md)
  integrates.

- [`rvt_insolation()`](https://wiesehahn.github.io/rvtr/reference/rvt_insolation.md):

  The **energy** rather than the duration, in kWh/m2, with beam and
  diffuse together. An hour of December sun is not an hour of June sun,
  so this is what
  [`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md)'s
  hour counts are really a proxy for.

They are complements rather than alternatives. Real ground receives
both, so a north-facing slope that
[`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md)
shows getting no sun at all is still lit diffusely, and a cell's total
light is neither number by itself. All three share the same horizon
search, so `reach`, `num_directions` and the multi-resolution levels of
[rvt_reach](https://wiesehahn.github.io/rvtr/reference/rvt_reach.md)
mean the same thing in each.

## See also

[`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md),
which integrates this over many sun positions instead of freezing one,
and
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md),
which shades by surface orientation without casting shadows at all.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_shadow(sun_elevation = 15, reach = 30)
```
