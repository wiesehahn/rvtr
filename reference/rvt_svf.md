# Sky-view factor

How much of the sky each cell can see, from 0 (completely enclosed - the
bottom of a narrow pit or gorge) to 1 (completely open - flat ground or
a summit). Flat ground gives exactly 1.

## Usage

``` r
rvt_svf(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  num_directions = 16,
  reach = 10,
  noise_removal = 0,
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
  mosaics are read across file boundaries, so tiles do not produce seams

- out_path:

  output GeoTIFF path (default: a temp file)

- num_directions:

  number of compass directions scanned (default 16)

- reach:

  look this far in **map units** (metres, normally) using a
  multi-resolution search, instead of `radius_max` pixels at full
  resolution. Cost then grows with the logarithm of the distance rather
  than in proportion to it, which is what makes a reach of hundreds of
  metres affordable on fine data. `radius_max` and `noise_removal` are
  ignored when this is set. See
  [rvt_reach](https://wiesehahn.github.io/rvtr/reference/rvt_reach.md).

- noise_removal:

  0-3; raises the minimum search radius to ignore near-cell noise (0,
  10, 20 or 40 percent of `reach`)

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

`out_path`, invisibly. `dem` is the first argument, so this is
pipe-friendly: `dem |> rvt_svf("svf.tif")`.

## Details

It looks like the terrain lit by a uniformly bright overcast sky. There
is no sun position, so nothing disappears into cast shadow and no
feature is favoured just because it happens to face the light - which
makes it a good general-purpose alternative to hillshade when features
may run in any direction. Hollows, ditches and the insides of banks come
out dark; ridges, mounds and open ground come out bright.

## How it works

From each cell the terrain is scanned outwards along `num_directions`
evenly spaced compass directions (16 by default, i.e. every 22.5
degrees):

- Along each direction, samples are taken out to `reach`, finely enough
  (three samples per pixel) that no pixel the ray crosses is skipped;
  duplicates are then dropped.

- For each sample, the angle from the cell up to that sample is computed
  from the elevation difference over the true ground distance. Pixel
  size is taken into account, so non-square pixels are handled
  correctly - but it also means **the elevation units must match the map
  units** (metres over metres) for the angles to be meaningful.

- The steepest angle found along a direction is that direction's
  horizon. Horizons below the horizontal are treated as horizontal,
  since a cell on a summit sees a full half-dome and no more.

- The result is the mean of `1 - sin(horizon angle)` over all
  directions.

Everything runs at the raster's native resolution, tile by tile, with
tiles overlapping by the search radius so the result never depends on
how the raster was split up.

## Tuning

- `reach` is the main control: it sets the size of the features the
  result responds to. It is in **map units** (metres, normally), so the
  same number means the same ground distance on any DEM. Short reaches
  emphasise fine local texture; long ones bring in shading from more
  distant landforms and take proportionally longer.

- `num_directions` trades smoothness for time (cost is roughly linear).
  16 is the usual choice; 32 reduces the faint directional banding that
  can appear in smooth terrain.

- `noise_removal` skips the innermost part of each ray. Raise it on
  noisy DEMs - vegetation or point-cloud speckle right next to a cell
  can otherwise set that cell's horizon on its own. Kokalj and Hesse
  note this markedly improves visibility of archaeological features
  where the data suffer from flight-strip misalignment or were gridded
  finer than the point density really supports.

## Recommended settings

Kokalj and Hesse (2017) suggest, for archaeological work on any terrain,
a **10 m** search radius with **16 directions** - which is the default
here, `reach = 10`, whatever the resolution of the DEM. Much larger
radii belong to other disciplines - they cite 10 km for meteorological
work - and only generalise archaeological detail away.

For display they suggest a linear stretch of 0.65-1.0 on varied terrain
and 0.9-1.0 on very flat ground, where the useful values are otherwise
squeezed into a narrow band just below 1.

It works best on moderate to steep topography, where it brings out
depressions and features lying on slopes. On flat or very gentle ground
it is largely limited to negative features - pits, ditches, quarries,
dolines - and becomes very sensitive to DEM noise, which is where
`noise_removal` earns its keep.

## References

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Prostor, kraj, čas 14.
Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

Zakšek, K., Oštir, K. and Kokalj, Ž. (2011) Sky-View Factor as a Relief
Visualization Technique. *Remote Sensing* 3, 398-415.

## See also

[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
for a related measure of the same horizon that spreads convex and
concave terrain over a wider range.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
rvt_svf(dem)
```
