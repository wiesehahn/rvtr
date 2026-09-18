# Multi-directional hillshade

One shaded-relief image lit from several directions at once, on a 0-1
scale. A single
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)
always loses whatever happens to run parallel to its light, and buries
one side of every feature in deep shadow. Lighting from several sides at
the same time keeps both: no direction is blind, and the dark side stays
readable.

## Usage

``` r
rvt_multi_hillshade(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  sun_elevation = 45,
  z_factor = 1,
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

- sun_elevation:

  height of the light sources above the horizon, in degrees, shared by
  all four (default 45)

- z_factor:

  vertical exaggeration applied to the elevations (default 1, no
  exaggeration)

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

The trade is contrast: averaging several lights fills the shadows in, so
the image is more even but flatter than a single hillshade, and
typically occupies a narrower range than 0-1. A percentile stretch when
plotting (`rvt_plot(p, minmax_pct_cut = c(2, 98))`) gets that contrast
back. Use the single-direction
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)
when you need a *known* light direction

- for a figure whose caption states it, or as the base layer of
  [`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md).

## How it works

A wrapper around GDAL's own `gdaldem hillshade -multidirectional` (via
[`gdalraster::dem_proc()`](https://firelab.github.io/gdalraster/reference/dem_proc.html)),
which implements the oblique-weighted method of Mark (1992): the surface
is lit from four sources, at azimuths 225, 270, 315 and 360 degrees, all
at `sun_elevation` above the horizon, and the four results are averaged
with weights that depend on the direction each cell faces - a source
lighting a slope square-on counts for more than one grazing it. The
single azimuth of
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)
therefore has no equivalent here; there is nothing to set.

GDAL computes this in one pass over 8-neighbour (Horn) derivatives, so
it is fast - roughly as quick as a single hillshade - and the result is
rescaled here from GDAL's 1-255 bytes to a single Float32 band in 0-1,
to match every other product in the package. Areas of NoData in the DEM
stay NoData.

## Tuning

Only `sun_elevation` and `z_factor` matter.

- `sun_elevation` behaves as in
  [`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md):
  lower picks out subtle relief more strongly, higher gives a flatter,
  more even image. 45 is GDAL's default and a good starting point - the
  weighting was designed around it - but 30-35 is worth trying when the
  result looks washed out. Below about 20 the image starts to look harsh
  even with four light sources.

- `z_factor` exaggerates the vertical scale. Raise it (2-5) for very
  flat ground; it is also the way to correct a DEM whose elevations are
  in different units from its coordinates.

## Recommended settings

Kokalj and Hesse (2017) recommend multi-directional hillshading over the
single-direction kind for general terrain interpretation, with the sun
elevation chosen for the terrain much as for a single hillshade: low on
flat ground, around 45 degrees in steep or complex topography.

Their own variant renders 16 separate directions and collapses them with
a principal components analysis. This does the same job in one band and
a fraction of the time; if you specifically need the PCA treatment, call
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)
once per azimuth and take those into your own PCA.

## References

Mark, R.K. (1992) *Multidirectional, oblique-weighted, shaded-relief
image of the Island of Dominica*. U.S. Geological Survey Open-File
Report 92-422.

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

## See also

[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)
for a single, named light direction.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_multi_hillshade()
```
