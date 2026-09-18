# Simple local relief model (SLRM)

Height of each cell above or below the smoothed terrain around it, in
the DEM's elevation units. Positive where the ground is locally raised
(banks, mounds, walls), negative where it is cut in (ditches, hollow
ways, pits), and 0 where it follows the local trend.

## Usage

``` r
rvt_slrm(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  radius = 20,
  tile_size = NULL,
  threads = rvt_threads(),
  overwrite = FALSE,
  progress = FALSE
)

rvt_tpi(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  radius = 20,
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

- radius:

  radius of the smoothing window, in **map units** (metres, normally;
  default 20)

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

This is the simplest of the trend-removal visualizations and still one
of the most used: subtracting a smoothed copy of the terrain from itself
strips out hillslopes and valleys, so low features read the same on a
slope as on level ground.
[`rvt_msrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_msrm.md)
does the same thing across a range of scales at once and is more
forgiving when you don't know the size of what you are looking for; this
one gives you a single, sharply defined scale that you set directly.

## How it works

The DEM is averaged over a square window of `radius`, and that smoothed
surface is subtracted from the original. What remains is everything
smaller than roughly the window - anything broader has been absorbed
into the average and cancels out.

The window average comes from a summed-area table, so cost does not grow
with the radius, and NoData cells are left out of the average rather
than counted as zero. Computation is at native resolution, tile by tile,
with tiles overlapping by `radius`.

## Tuning

- `radius` is the one control, and it is in **map units**: the default
  20 is 20 m on a 1 m DEM but only 5 m on a 0.25 m DEM, so scale it with
  your resolution. It should be comfortably larger than the features you
  want to keep - features approaching the window size get partly
  absorbed into the trend and fade out.

- Too small a radius flattens everything and leaves only noise; too
  large leaves the hillslopes in and defeats the point. Somewhere around
  two to three times the width of your target features is a reasonable
  starting point.

- A square averaging window leaves faint blocky artefacts around very
  sharp features. If that becomes distracting,
  [`rvt_msrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_msrm.md)
  averages several scales and largely avoids it.

## Recommended settings

Kokalj and Hesse (2017) suggest a filter radius of about **10 m**
generally, 5 m on very flat ground and 25 m in steep or complex
terrain - metres again, and `radius` is in map units, so the same number
means the same distance whatever the resolution. For display they use a
linear stretch of -1 to +1 m, tightening to -0.5 to +0.5 m on very flat
ground and opening to -2 to +2 m on steep ground; the right choice
really depends on how tall the features you are after are.

Trend removal is picked out as one of the most useful approaches for
very low relief on flat to moderate ground - former field boundaries,
levelled burial mounds. Compared with the closely related local relief
model, this simple form is faster and simpler, at the cost of somewhat
less faithful relative heights for the anomalies it reveals.

## Also known as topographic position index

With a square window, "cell minus the mean around it" is exactly Weiss's
**Topographic Position Index**, so `rvt_tpi()` is provided as an alias
and the two compute the same thing. The names come from different
literatures - TPI from landform classification, SLRM from archaeological
prospection - and people search for one or the other.

`gdalraster::dem_proc(mode = "TPI")` also computes TPI but is fixed to a
3x3 window; this takes any `radius`, at the same cost thanks to the
summed-area table (`radius = 1` reproduces the 3x3 case). For the
standardised version, dividing by the local standard deviation, see
[`rvt_dev()`](https://wiesehahn.github.io/rvtr/reference/rvt_dev.md).

## References

Hesse, R. (2010) LiDAR-derived Local Relief Models - a New Tool for
Archaeological Prospection. *Archaeological Prospection* 17, 67-72.

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

## See also

[`rvt_msrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_msrm.md)
for the multi-scale version,
[`rvt_dev()`](https://wiesehahn.github.io/rvtr/reference/rvt_dev.md) for
the standardised one.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_slrm()
```
