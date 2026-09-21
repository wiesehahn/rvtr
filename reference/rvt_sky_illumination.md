# Sky illumination

How much diffuse daylight each cell receives from the sky as a whole -
accounting both for the way the ground is tilted and for the surrounding
terrain blocking part of the view. Scaled so flat, open ground reads 1;
below that means shaded by slope or surroundings.

## Usage

``` r
rvt_sky_illumination(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  sky_model = c("overcast", "uniform"),
  num_directions = 32,
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

- sky_model:

  `"overcast"` (default) or `"uniform"` - how sky brightness varies with
  height above the horizon

- num_directions:

  number of directions searched (default 32)

- reach:

  how far to look, in **map units** (metres, normally; default 100).
  Scanned at full resolution for the first `pyramid_px` cells (100 by
  default, so 100 m on a 1 m DEM); past that the search reads from
  progressively coarser copies of the DEM, so a long reach stays
  affordable - see
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

`out_path`, invisibly

## Details

Unlike
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)
there is no sun and no cast shadow, and unlike
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) the
sky is not counted uniformly: what a surface actually receives depends
on how squarely it faces each part of the sky, so a slope tilted away
from open sky is darker even with an identical horizon. The result looks
like an overcast day, which is a fair approximation of how terrain reads
to the eye in soft light.

This is the most expensive metric here by a wide margin - see Tuning.

## How it works

The horizon is found exactly as in
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) - a
maximum slope along each of `num_directions` rays out to `reach`. Each
direction then contributes light according to how much sky it leaves
open above the horizon *and* how that sky sits relative to the surface:
the cell's slope and aspect weight every direction, so sky behind the
surface contributes nothing. The `sky_model` sets how brightness varies
with height in the sky

- `"uniform"` treats it as equally bright everywhere, `"overcast"` makes
  it brighter overhead, which is the standard cloudy-sky distribution
  and the default.

Results are scaled by the value flat ground produces, so flat reads 1
regardless of `sky_model` or `num_directions`.

## Tuning

- `reach` dominates the cost: the default 100 m means around 3500
  terrain samples per cell, some 20 times what
  [`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md)
  does at its own default. Expect around a minute for a 4000 x 4000
  raster. Shortening it to 30 m still captures the terrain that matters
  most for local shading; going *longer* is much cheaper than it used to
  be, since past one full-resolution scan - `pyramid_px` cells, 100 by
  default - the search coarsens as it goes (see
  [rvt_reach](https://wiesehahn.github.io/rvtr/reference/rvt_reach.md)).

- `num_directions` costs time in direct proportion. 32 is the default;
  16 halves the work at some loss of smoothness.

- `sky_model = "uniform"` is slightly cheaper and gives a flatter, more
  even image; `"overcast"` gives more contrast between open and enclosed
  ground.

## Recommended settings

Kokalj and Hesse (2017) mention diffuse insolation only in passing, and
their verdict is worth repeating before you spend the runtime: it does
reveal archaeological features, but it needs considerably more
computation than the alternatives and the results come out **more
generalised** than, say, sky-view factor. If you want a diffuse-light
view of the terrain and are not specifically after the physical
illumination quantity,
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md)
gets you most of the way for a twentieth of the work, and
[`rvt_asvf()`](https://wiesehahn.github.io/rvtr/reference/rvt_asvf.md)
adds back a sense of direction.

## Direct sun versus diffuse sky

Three functions here measure light reaching the ground, and they answer
different questions:

- [`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md):

  **Direct** sunlight, **integrated over time**: how many hours the sun
  itself reaches a cell across a day, a season or a year. Needs a date
  range and a latitude, and answers in hours. A clear-sky quantity.

- `rvt_sky_illumination()`:

  **Diffuse** skylight, at **no particular time**: how much light a cell
  receives from the sky dome as a whole, given both its horizon and the
  way it tilts. No date, no sun position, no latitude; the answer is
  relative, with flat open ground at 1. An overcast-sky quantity.

- [`rvt_shadow()`](https://wiesehahn.github.io/rvtr/reference/rvt_shadow.md):

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

## References

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

## See also

[`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md)
for direct sun over time rather than diffuse light at no particular
time, and
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) for
a much cheaper measure of the same horizon that ignores which way the
ground faces.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
# small radius to keep the example quick
dem |> rvt_sky_illumination(reach = 10, num_directions = 8)
```
