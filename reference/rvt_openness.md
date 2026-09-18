# Positive openness

The average angle from straight up (the zenith) down to the horizon, in
degrees. Flat ground gives exactly 90 degrees. Above 90 means the cell
stands proud of its surroundings - ridges, mounds, bank tops, wall
lines; below 90 means it sits enclosed - ditches, hollows, the inside of
a bank.

## Usage

``` r
rvt_openness(
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
  mosaics are read across file boundaries, so tiles do not produce
  seams. A terra `SpatRaster` works too: one read from a file is used as
  it lies, while one terra computed in memory is written to a temporary
  GeoTIFF first.

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

This measures the same horizon as
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) but
keeps more of its range. Sky-view factor compresses everything above the
horizontal into a narrow band just below 1, whereas openness spreads
convex, flat and concave terrain out evenly around 90 degrees - so the
edges of low earthworks tend to come out crisper. Values are not
restricted to 0-90: around 107 degrees is reached on strongly convex
terrain in the bundled sample tile.

## How it works

The horizon search is identical to
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) -
same directions, same sampling along each ray, same use of true ground
distance - with two differences in how the horizon angles are combined:

- Horizons below the horizontal are **not** clamped. A cell looking down
  into falling ground keeps its negative angle, which is what lets
  convex terrain rise above 90 degrees.

- The result is `90 - mean(horizon angle)` in degrees, over all
  directions.

Like sky-view factor it runs at native resolution, tile by tile, with
tiles overlapping by the search radius. Elevation units must match the
map units for the angles to mean anything.

## Tuning

The parameters behave exactly as in
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) -
`reach` sets the size of features the result responds to, in map units,
`num_directions` trades smoothness against time, and `noise_removal`
protects against speckle immediately next to a cell.

A common pairing is to compute this alongside
[`rvt_openness_negative()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness_negative.md):
the positive image favours convex features, the negative one concave,
and comparing or differencing the two separates banks from ditches more
clearly than either does alone.

## Recommended settings

As for
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md):
Kokalj and Hesse (2017) use a **10 m** radius (`reach = 10`) with 16
directions, and suggest displaying positive openness with a linear
stretch of 65-95 degrees on varied terrain, 85-91 on very flat ground,
55-95 on steep or complex ground.

Because openness ignores the mathematical horizon and takes in the whole
sphere, the image comes out notably *flat* - the broad shape of the
landscape largely disappears, rather as if the trend had been removed.
That costs you the intuitive read of the terrain that sky-view factor
keeps, and makes interpretation a little harder, but it buys two things:
the result is not washed out by gentle or steep slopes and so works
across varied topography, and a given feature looks much the same
wherever it sits

- on the flat or partway up a hillside. Doneus (2013) makes the point
  that this consistency of "signature" is what makes openness well
  suited to automated feature detection, not just to looking.

## References

Yokoyama, R., Shirasawa, M. and Pike, R. J. (2002) Visualizing
Topography by Openness: A New Application of Image Processing to Digital
Elevation Models. *Photogrammetric Engineering and Remote Sensing* 68,
251-266.

Doneus, M. (2013) Openness as Visualization Technique for Interpretative
Mapping of Airborne LiDAR Derived Digital Terrain Models. *Remote
Sensing* 5, 6427-6442.

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

## See also

[`rvt_openness_negative()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness_negative.md)
for the concave counterpart,
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) for
the 0-1 sky-fraction version of the same horizon.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_openness()
```
