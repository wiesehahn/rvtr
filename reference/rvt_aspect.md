# Aspect

The compass direction each cell faces - the way water would run off it.
0 is north, 90 east, 180 south, 270 west.

## Usage

``` r
rvt_aspect(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  units = c("degree", "radian"),
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

  `"degree"` (default) or `"radian"`

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

Aspect on its own is rarely a useful picture: it wraps around at north,
so a smooth north-facing slope shows a hard seam between 359 and 0, and
flat ground points in an essentially arbitrary direction. It earns its
place as an input to other things - solar or ecological modelling, or
splitting a slope map into the sides of a feature.

## How it works

Taken from the same finite-difference derivatives as
[`rvt_slope()`](https://wiesehahn.github.io/rvtr/reference/rvt_slope.md),
as the azimuth of steepest descent, then wrapped into 0-360. Cells on
perfectly flat ground have no meaningful aspect; the underlying formula
carries a small numerical guard so they come out as a fixed direction
rather than undefined, which is worth knowing before you read anything
into large uniform patches.

## Tuning

Nothing to tune. If you want a *readable* directional image rather than
raw azimuths,
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)
or
[`rvt_asvf()`](https://wiesehahn.github.io/rvtr/reference/rvt_asvf.md)
are what you actually want - both fold aspect into a shaded result
without the wrap-around discontinuity.

Note this will not match `gdalraster::dem_proc(mode = "aspect")`
exactly. That uses Horn's eight-neighbour estimate; this uses the
four-neighbour central difference, so that slope and aspect here stay
consistent with each other and with
[`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md). On
the bundled tile the two agree to about 1.4 degrees on slopes above 20
degrees but diverge to 33 degrees on ground below 1 degree - which is
not a disagreement about the terrain so much as aspect being ill-defined
when there is barely a gradient to point along.

## See also

[`rvt_slope()`](https://wiesehahn.github.io/rvtr/reference/rvt_slope.md),
which shares the same derivative pass.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
rvt_aspect(dem)
```
