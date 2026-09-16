# Anisotropic sky-view factor

Sky-view factor computed under a sky that is brighter on one side,
instead of uniformly bright all round. Still 0-1, and still 1 on flat
ground, but terrain now casts a soft directional shading: slopes facing
the bright side come out lighter, slopes turned away from it darker.

## Usage

``` r
rvt_asvf(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  num_directions = 16,
  reach = 10,
  noise_removal = 0,
  main_direction = 315,
  level = 1,
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

- main_direction:

  compass azimuth, in degrees, of the brightest part of the sky (default
  315, i.e. north-west)

- level:

  strength of the anisotropy: 1 gentle (default), 2 strong

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

The point is to recover some of the readability of a hillshade - the eye
is good at reading relief lit from one side - without a hillshade's main
flaw, which is that features running parallel to the light almost
disappear. Nothing is ever fully in shadow here, so a feature aligned
with `main_direction` is dimmed rather than lost.

## How it works

The horizon search is exactly the one in
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md).
The only difference is in the final averaging: instead of every
direction counting equally, each gets a weight running from 1 in
`main_direction` down to `min_weight` at the opposite azimuth, following
`cos(half the angle from main_direction)` raised to a power. The
weighted contributions are divided by the total weight, which keeps the
result on the same 0-1 scale as plain sky-view factor.

`level` picks how sharply the weight falls away: level 1 is gentle
(opposite sky still counts 0.4), level 2 is strong (opposite sky counts
0.1), which gives a more hillshade-like image.

## Tuning

- `main_direction` is a compass azimuth in degrees - 0 north, 90 east,
  180 south, 270 west. The 315 default is north-west, the conventional
  illumination direction for terrain images, and matches
  [`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)'s
  default so the two can be compared or blended. Rotate it to bring out
  features running in a particular direction: relief perpendicular to
  the light shows up most.

- `level` 1 keeps the result close to plain sky-view factor with a hint
  of direction; 2 pushes the directional effect much harder.

- `reach`, `num_directions` and `noise_removal` behave exactly as in
  [`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md).
  Raising `num_directions` matters a little more here than for plain
  sky-view factor, since the weights vary between directions and a
  coarse set samples that variation coarsely.

## Recommended settings

Radius and directions as for
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) -
Kokalj and Hesse (2017) use 10 m and 16 directions for archaeological
work, and the same 0.65-1.0 / 0.9-1.0 display stretches apply.

Their assessment is that the anisotropy "brings back some of the
plasticity of hill shading and gives better details on very flat areas",
which is a fair summary of when to prefer this over plain sky-view
factor: reach for it when the isotropic image reads as flat and hard to
interpret, particularly on level ground. Both are described as very good
general- purpose choices, because they show small features whatever
their orientation and shape, on most terrain.

## References

Zakšek, K., Oštir, K., Pehani, P., Kokalj, Ž. and Polert, E. (2012) Hill
Shading Based on Anisotropic Diffuse Illumination. In *Symposium GIS
Ostrava 2012*, 1-10. Ostrava: Technical University of Ostrava.

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

## See also

[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) for
the uniform-sky version,
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)
for a conventional directional shading.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
rvt_asvf(dem)
```
