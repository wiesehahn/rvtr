# Geomorphons

Classifies every cell into one of ten landform types - flat, peak,
ridge, shoulder, spur, slope, hollow, footslope, valley, pit - by
looking at the shape of the terrain around it rather than by measuring a
single quantity.

## Usage

``` r
rvt_geomorphons(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  reach = 20,
  skip = 1,
  flat_threshold = 1,
  flat_distance = 0,
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

- reach:

  how far to look, in **map units** (metres, normally; default 20)

- skip:

  innermost distance to ignore, in map units; 0 skips nothing (default
  1)

- flat_threshold:

  angle below which ground counts as flat, in degrees (default 1)

- flat_distance:

  distance beyond which the flatness threshold is relaxed, in map units;
  0 disables it (default 0)

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

`out_path`, invisibly - class codes, see
[rvt_geomorphon_classes](https://wiesehahn.github.io/rvtr/reference/rvt_geomorphon_classes.md)

## Details

This is the odd one out in the package: the output is **categorical**,
not a continuous surface. Values are the class codes in
[rvt_geomorphon_classes](https://wiesehahn.github.io/rvtr/reference/rvt_geomorphon_classes.md),
and it should be displayed with a discrete palette, never a continuous
ramp - the numbers are labels, and "6" is not between "5" and "7" in any
meaningful sense.

What makes it useful is that the classification is *scale-explicit and
self-normalising*: because a cell is judged by the angles to what it can
see rather than by absolute heights, the same landform gets the same
label whether it sits in a mountain range or on a floodplain. That makes
it a good starting point for asking "where are all the ridges" or "where
are all the hollows" across ground of mixed character, and a common
input to further analysis rather than a picture in its own right.

## How it works

Along each of eight directions the algorithm looks out to `reach` and
records the highest and lowest line-of-sight angles - the most it has to
look up, and the most it has to look down. Whichever is larger decides
that direction's verdict: terrain rises away, falls away, or is level
within the `flat_threshold`. That produces eight ternary votes.

Counting how many directions rise and how many fall gives a pair of
numbers that indexes a lookup table of forms. Eight falling directions
is a peak, eight rising is a pit, half and half is a slope, and the rest
of the table fills in ridges, spurs, hollows, footslopes and so on. Only
the counts matter, not which particular directions - so the
classification is rotation-invariant by construction.

It is the same ray walk as
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md), so
cost scales the same way with `reach`; it just keeps both extremes per
direction instead of the maximum alone. Those two extremes are exactly
what
[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
and
[`rvt_openness_negative()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness_negative.md)
average - geomorphons was introduced as an extension of terrain
openness, classifying the pattern of those angles rather than averaging
them away.

## Tuning

- `reach` is the main control and sets **the scale of landform you
  get**. Small values classify micro-relief - individual banks and
  ditches; large values classify the hill the whole site sits on. It is
  in **map units**, so scale it with your resolution. There is no
  universally right value: choose it from the size of the landforms you
  care about.

- `flat_threshold` (degrees) decides how level counts as flat. Raise it
  to sweep gently undulating ground into the flat class, lower it to
  sub-divide subtle terrain. On noisy DEMs too small a value produces
  speckle.

- `skip` ignores the innermost ring, exactly as `noise_removal` does for
  [`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md),
  and is the right first response to a speckled result.

- `flat_distance` (map units, 0 to disable) relaxes the flatness
  threshold with range: past this distance a fixed angle would make
  every distant slope significant, so beyond it the threshold becomes
  the angle the same height difference subtends. Worth setting when
  `reach` is large.

## References

Jasiewicz, J. and Stepinski, T. (2013) Geomorphons - a pattern
recognition approach to classification and mapping of landforms.
*Geomorphology* 182, 147-156.
[doi:10.1016/j.geomorph.2012.11.005](https://doi.org/10.1016/j.geomorph.2012.11.005)

## See also

[rvt_geomorphon_classes](https://wiesehahn.github.io/rvtr/reference/rvt_geomorphon_classes.md)
for the class codes;
[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
and
[`rvt_openness_negative()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness_negative.md)
for the continuous measures this classifies.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_geomorphons(reach = 10)
```
