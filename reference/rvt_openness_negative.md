# Negative openness

The same measurement as
[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md),
made on the terrain turned upside down. Concave features become the
prominent ones: ditches, pits, hollow ways and quarry scoops read high,
where positive openness favours banks, mounds and ridges.

## Usage

``` r
rvt_openness_negative(
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

`out_path`, invisibly

## Details

Note that only the *shape* preference is inverted, not the numbers. Flat
ground still gives exactly 90 degrees, and a cell in a pit reads *above*
90 here - because on the inverted surface that pit is a mound. Values
are in degrees and, as with positive openness, are not restricted to
0-90.

## How it works

Elevations are negated, and then exactly the horizon search described in
[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
runs on the result - same directions, same sampling, same
`90 - mean(horizon angle)`. Nothing else differs, so the two are
directly comparable and can be differenced.

## Tuning

Identical to
[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md).
Use the same `reach` for both if you intend to compare or difference
them - a mismatch there makes the two images respond to different
feature sizes and the comparison stops meaning much.

## Recommended settings

Worth being clear about, since the name invites the assumption:
**negative openness is not the inverse of positive openness.** It
carries genuinely different information rather than the same image with
the sign flipped. Where positive openness picks out the convexities -
the ridge *between* two hollow ways, the raised rim of a crater -
negative openness picks out the deepest parts of the concavities: the
hollow way itself, the floor of a gorge, the foot of a cliff. Computing
both and reading them together is the point.

Kokalj and Hesse (2017) recommend displaying it with an **inverted**
greyscale, dark for high values, so that concave features read dark in
both images and the pair stay visually consistent. Their suggested
stretch is 60-95 degrees on varied terrain, 75-95 on very flat ground
and 45-95 on steep or complex ground. Radius and directions as for
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md): 10
m, 16 directions.

## References

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

## See also

[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
for the convex counterpart.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
rvt_openness_negative(dem)
```
