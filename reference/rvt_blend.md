# Blend a layer over a raster

Lays one raster over another with a blend mode, building up a display
image the way layers stack in QGIS. Pipe rasters through it and render
once at the end.

## Usage

``` r
rvt_blend(x, layer, mode = "normal", opacity = 1, range = NULL, invert = FALSE)
```

## Arguments

- x:

  a raster path, or an `rvt_stack` to add to

- layer:

  path to the raster to lay over `x`

- mode:

  one of
  [rvt_blend_modes](https://wiesehahn.github.io/rvtr/reference/rvt_blend_modes.md);
  `"normal"` (default) just replaces, which with `opacity` is a plain
  fade between the two

- opacity:

  0-1, applied after the mode (default 1)

- range:

  `c(lo, hi)` to stretch `layer` over, or `NULL` (default) for its own
  minimum and maximum

- invert:

  flip the layer after stretching, so high reads dark. Slope and
  negative openness are usually shown this way.

## Value

an `rvt_stack`

## Details

    hillshade |>
      rvt_blend(openness, "overlay",  opacity = 0.5) |>
      rvt_blend(svf,      "multiply", opacity = 0.25) |>
      rvt_render("relief.tif")

**Order is bottom to top.** The raster piped in is the *background*; the
one passed as an argument goes *over* it. That matters, because
`overlay`, `dodge`, `burn`, `hard_light` and `subtract` all treat their
two inputs differently. It reads the way a QGIS layer panel does,
upwards.

**Nothing is computed until you render.** `rvt_blend()` returns a
recipe, not a raster, so a chain of layers costs one pass over the data
instead of one per layer.
[`rvt_render()`](https://wiesehahn.github.io/rvtr/reference/rvt_render.md)
writes the file;
[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)
will draw a stack directly, so exploring needs no render call at all.

## Choosing a mode

Three of the thirteen modes cover most composites. `multiply` only ever
darkens, so it is what you want for deepening enclosed ground with a
sky-view factor. `screen` only ever lightens. `overlay` does both at
once - darkening what is already dark and lightening what is already
light - which adds contrast without washing the image out, and is why
[`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md)
uses it for its openness layer.

The rest are situational: `darken`/`lighten` take the extreme of the
two, `difference` shows where two visualizations disagree, `soft_light`
is a gentler `overlay`, and `dodge`/`burn` are aggressive enough to
clip. All thirteen are in
[rvt_blend_modes](https://wiesehahn.github.io/rvtr/reference/rvt_blend_modes.md).

## Ranges

Every layer is stretched to 0-1 before blending. Leave `range` alone and
the raster's own minimum and maximum are used, measured once before any
tile is read. Set it to concentrate contrast where it matters -
`c(0.7, 1)` on a sky-view factor spreads out the band where ordinary
terrain actually lives.
[`rvt_range()`](https://wiesehahn.github.io/rvtr/reference/rvt_range.md)
reads a percentile stretch off the raster.

## Different resolutions

Layers may differ in resolution by a whole factor - a 0.2 m orthophoto
under a 1 m hillshade, say - as long as they cover the same extent in
the same coordinate system. The result is rendered at the finest
resolution in the stack, and coarser layers are upsampled onto it with
bilinear interpolation when rendering. Rasters fetched for one place
with
[`rvt_data_lgln()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln.md)
or
[`rvt_data_mapterhorn()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_mapterhorn.md)
line up this way. Anything else stops with an error naming the mismatch,
rather than blending misaligned.

A three-band layer is stretched over a single range spanning all its
bands, never one range per band: stretching each channel to its own
extremes would shift the colour balance. Stretch the channels yourself
first if that is what you want.

## See also

[`rvt_render()`](https://wiesehahn.github.io/rvtr/reference/rvt_render.md),
[`rvt_range()`](https://wiesehahn.github.io/rvtr/reference/rvt_range.md),
[rvt_blend_modes](https://wiesehahn.github.io/rvtr/reference/rvt_blend_modes.md);
[`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md) for
a ready-made archaeological composite.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
hs <- dem |> rvt_hillshade()
svf <- dem |> rvt_svf()
hs |> rvt_blend(svf, "multiply", opacity = 0.25)
#> <rvt_stack> 2 layers, bottom to top
#>   1. file223e11bb4f14.tif         (background) opacity 1.00  range auto
#>   2. file223e6edf8135.tif         multiply    opacity 0.25  range auto
#> Renders at 1000 x 1000, 1 m. Not a raster yet - call rvt_render() to write it.
```
