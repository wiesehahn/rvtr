# Hillshade

The familiar shaded-relief image: terrain lit by a single low sun,
bright where slopes face the light and dark where they turn away, from 0
to 1.

## Usage

``` r
rvt_hillshade(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  sun_azimuth = 315,
  sun_elevation = 35,
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

- sun_azimuth:

  compass azimuth of the sun, in degrees (default 315, i.e. north-west)

- sun_elevation:

  height of the sun above the horizon, in degrees (default 35)

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

It is the most intuitive way to look at terrain, which is why it remains
the default in most GIS software - but it has one well-known weakness
worth keeping in mind: features running parallel to the light direction
almost vanish, while those across it are exaggerated. If that matters,
either render a few different `sun_azimuth` values and compare, or use a
direction-free measure such as
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md), or
[`rvt_asvf()`](https://wiesehahn.github.io/rvtr/reference/rvt_asvf.md)
which keeps a sense of direction without letting anything disappear
entirely.

## How it works

The slope and the direction each cell faces are estimated from its
immediate neighbours (as in
[`rvt_slope()`](https://wiesehahn.github.io/rvtr/reference/rvt_slope.md)),
then combined with the sun position into the cosine of the angle between
the surface and the incoming light. Values that come out negative -
surfaces pointing away from the sun

- are clamped to 0.

This is local shading only: it says whether a surface is tilted towards
the light, not whether anything stands between it and the sun, so no
long shadows are cast across the terrain. Computation is at native
resolution, tile by tile.

## Tuning

- `sun_azimuth` is a compass azimuth in degrees - 0 north, 90 east, 180
  south, 270 west. The 315 (north-west) default is a long-standing
  convention: lighting from the upper left is what makes relief read as
  raised rather than sunken to most viewers. Rotating it to run across
  the grain of the features you care about brings them out best.

- `sun_elevation` is the sun's height above the horizon in degrees.
  Lower values (15-25) exaggerate faint relief and are the usual choice
  for subtle archaeological features; higher values give a flatter, more
  even image.

## Recommended settings

Kokalj and Hesse (2017) tie the sun elevation to the terrain: about **35
degrees** generally, **below 10** to draw out low relief on flat or
gently sloping ground, and **above 45** in steep or complex topography,
where a low sun simply burns out one side of every hill and blacks out
the other. Pushed towards vertical on steep ground the image starts to
resemble
[`rvt_slope()`](https://wiesehahn.github.io/rvtr/reference/rvt_slope.md),
which they suggest as the better tool there. Azimuth stays at 315
throughout; display with a linear stretch and a 2% cut-off.

They recommend starting any investigation with a shaded relief overview,
because it gives the most natural-looking read of the topography and
helps you judge which other techniques are worth applying - but also
warn about three failure modes: features running parallel to the light
nearly vanish (so vary the azimuth, or use
[`rvt_multi_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_multi_hillshade.md)),
slopes facing into or away from the light saturate to white or black,
and relief can read *inverted*, hollows appearing as mounds. That last
illusion is why lighting from the north-west is conventional: most
viewers read top-left lighting correctly, and other azimuths invite the
illusion.

## References

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

## See also

[`rvt_asvf()`](https://wiesehahn.github.io/rvtr/reference/rvt_asvf.md)
for a directional image where nothing is lost in shadow,
[`rvt_slope()`](https://wiesehahn.github.io/rvtr/reference/rvt_slope.md)
for the underlying steepness.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_hillshade()
```
