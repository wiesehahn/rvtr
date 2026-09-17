# Save a raster or blend as an image

Writes a WebP or JPEG picture of any raster or blend stack, at full
resolution: a colour palette and stretch for a single band, the colours
as they are for three bands. Where
[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)
draws a quick look on screen, this makes the file you put on a web page
or in a report.

## Usage

``` r
rvt_image(
  x,
  out_path,
  col = "Grays",
  range = NULL,
  band = NULL,
  quality = 90,
  tile_size = NULL,
  threads = rvt_threads(),
  overwrite = FALSE,
  progress = FALSE
)
```

## Arguments

- x:

  a raster path, or an `rvt_stack` from
  [`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)

- out_path:

  the image to write, ending in `.webp` or `.jpg`

- col:

  palette name or colours for a single band, as in
  [`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)
  (default `"Grays"`); see
  [`rvt_palettes()`](https://wiesehahn.github.io/rvtr/reference/rvt_palettes.md)

- range:

  `c(lo, hi)` to stretch a single-band or RGB raster over, or `NULL`
  (default) for its own minimum and maximum

- band:

  which band to draw, or three band numbers for red, green and blue.
  Defaults to band 1, or 1-3 for a three-band raster.

- quality:

  compression quality, 1-100 (default 90)

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

The format follows the file extension: `.webp` or `.jpg`. Images carry
no georeferencing; for a result to use in a GIS, keep the `.tif` the
metric function wrote.

[`rvt_render()`](https://wiesehahn.github.io/rvtr/reference/rvt_render.md),
[`rvt_relight()`](https://wiesehahn.github.io/rvtr/reference/rvt_relight.md)
and [`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md)
write images directly too - give them an image path - so their results
need no second step.

## Choosing a format

Both are compressed for viewing, not for keeping values: the default
quality of 90 looks the same as the original, but the numbers are not
preserved. Keep the GeoTIFF for anything that needs them.

**WebP** is the smaller file in every case measured - about 10% smaller
than JPEG for photographs and shading, and a tenth of a compressed
GeoTIFF - and shows NoData as transparent.

**JPEG** is the one to reach for when a tool cannot read WebP, or for
images wider than WebP's limit of 16,383 pixels (JPEG allows 65,500). It
has no transparency, so NoData is black.

## Stretch and colours

A single band is stretched across `range` and coloured with `col`,
exactly as in
[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md).
`range = NULL` uses the raster's own minimum and maximum;
[`rvt_range()`](https://wiesehahn.github.io/rvtr/reference/rvt_range.md)
gives a percentile stretch, and each metric's help page lists
recommended ranges. The range is fixed before any tile is read, so the
picture does not depend on `tile_size`.

Three bands - an RGB result such as
[`rvt_mstp()`](https://wiesehahn.github.io/rvtr/reference/rvt_mstp.md),
or a blend stack over an orthophoto - are stretched over one range
spanning all bands, and `col` is not used. A blend stack is stretched
layer by layer when it is blended (see
[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)),
so `range` does not apply to it either.

## See also

[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)
for a quick look on screen,
[`rvt_render()`](https://wiesehahn.github.io/rvtr/reference/rvt_render.md)
for blends

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
svf <- rvt_svf(dem)
rvt_image(svf, tempfile(fileext = ".webp"), range = c(0.7, 1))

# signed data: a diverging palette with a symmetric range
dem |> rvt_slrm() |>
  rvt_image(tempfile(fileext = ".webp"), "Blue-Red 3", range = c(-1, 1))
```
