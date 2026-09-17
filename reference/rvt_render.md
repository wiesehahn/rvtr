# Render a blend stack to a raster

Writes the composed image built by
[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)
as a Cloud-Optimized GeoTIFF, values 0-1, ready to display - or straight
to a WebP or JPEG picture when `out_path` ends in `.webp` or `.jpg`.

## Usage

``` r
rvt_render(
  stack,
  out_path = fs::file_temp(ext = "tif"),
  quality = 90,
  tile_size = NULL,
  threads = rvt_threads(),
  overwrite = FALSE,
  progress = FALSE
)
```

## Arguments

- stack:

  an `rvt_stack` from
  [`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)

- out_path:

  where to write; defaults to a temporary GeoTIFF

- quality:

  for `.webp` and `.jpg`, compression quality 1-100 (default 90)

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

This is where the work happens. Layers are read tile by tile and blended
in one pass, so a stack of any depth costs a single traverse of the
data. Each layer's stretch range is measured once, before the first
tile, which is what keeps the result independent of `tile_size`.

## Output

The file extension decides the format. `.tif` writes a Cloud-Optimized
GeoTIFF to use in a GIS. `.webp` and `.jpg` write a plain picture with
no georeferencing, for web pages and reports, with no GeoTIFF made along
the way; see
[`rvt_image()`](https://wiesehahn.github.io/rvtr/reference/rvt_image.md)
for how the two formats compare. A single-band blend is written in grey;
for a colour palette use `rvt_image(stack, "out.webp", col = ...)`.

## See also

[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
rvt_hillshade(dem) |>
  rvt_blend(rvt_svf(dem), "multiply", opacity = 0.25) |>
  rvt_render()
```
