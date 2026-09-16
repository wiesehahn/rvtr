# Slope

Steepness of the terrain at each cell, ignoring which way it faces. Flat
ground is 0 and a vertical face is 90 degrees (or 1.571 radians, or an
unbounded percentage).

## Usage

``` r
rvt_slope(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  units = c("degree", "radian", "percent"),
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

- units:

  `"degree"` (default), `"radian"`, or `"percent"` (100 times the rise
  over the run, so 45 degrees is 100 percent)

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

On its own this is a plain terrain derivative rather than a relief
visualization, but it is a useful layer to have: it picks out breaks of
slope very sharply, so the edges of banks, terraces, ditches and quarry
faces stand out as bright lines even where the features themselves are
low. It is also one of the ingredients
[`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md)
blends.

## How it works

Slope comes from a central difference over the **four** orthogonal
neighbours: the east-west gradient is taken from the two horizontal ones
and the north-south gradient from the two vertical ones, each divided by
the corresponding pixel size, and the slope is the arctangent of their
combined magnitude.

Several estimators are in common use and they do not give the same
answer, so it is worth naming the one here. This four-neighbour form
responds to the finest detail in the DEM, and to the most noise; Horn's
eight-neighbour estimate (what `gdaldem` and
[`gdalraster::dem_proc()`](https://firelab.github.io/gdalraster/reference/dem_proc.html)
use) and the Evans-Young quadratic fit average over a wider stencil and
come out smoother. The four-neighbour form is used here so that slope,
aspect and hillshade stay consistent with each other and with
[`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md).

Because horizontal distances come from the raster's own resolution, the
elevation units must match the map units (metres over metres) for the
angles to be right. Computation is at native resolution, tile by tile.

## Tuning

There is nothing to tune but the output units. Slope has no radius or
smoothing parameter: it is always measured across single pixels, which
means a noisy DEM gives a noisy slope. If that is a problem, smooth the
DEM before calling this, or reach for a metric with a built-in scale
such as
[`rvt_msrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_msrm.md)
or
[`rvt_local_dominance()`](https://wiesehahn.github.io/rvtr/reference/rvt_local_dominance.md).

## Recommended settings

Kokalj and Hesse (2017) recommend displaying slope in **inverted**
greyscale - steep dark, flat white - which keeps a surprisingly plastic,
relief-like impression of the terrain. Their suggested linear stretches
are 0-50 degrees generally, 0-15 on very flat ground, 0-60 on steep or
complex ground. They note Challis et al. found slope the single best
technique across most of the situations they tested.

Its one serious limitation is worth stating plainly: **slope cannot tell
a bank from a ditch.** Rising and falling ground of the same gradient
are drawn identically, so a profile across a feature can easily be read
as convex when it is really concave. Pair it with something that carries
sign

- [`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
  and
  [`rvt_openness_negative()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness_negative.md),
  or
  [`rvt_msrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_msrm.md),
  whose values are signed by construction - before committing to an
  interpretation.

## References

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

Pike, R. J., Evans, I. S. and Hengl, T. (2009) Geomorphometry: a brief
guide. In: Hengl, T. and Reuter, H. I. (eds.) *Geomorphometry: Concepts,
Software, Applications*. Developments in Soil Science 33. Amsterdam:
Elsevier, 3-30.
[doi:10.1016/S0166-2481(08)00001-9](https://doi.org/10.1016/S0166-2481%2808%2900001-9)

## See also

[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md),
which combines slope with the direction the terrain faces.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
rvt_slope(dem)
```
