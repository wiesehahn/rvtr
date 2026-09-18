# Start a blend stack

The background layer a stack is built on. Only needed when the
background itself needs a `range` or `invert`; otherwise pipe the raster
straight into
[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md),
which starts a stack for you.

## Usage

``` r
rvt_stack(base, range = NULL, invert = FALSE)
```

## Arguments

- base:

  path to the background raster

- range:

  `c(lo, hi)` to stretch over, or `NULL` (default) for the raster's own
  minimum and maximum. See
  [`rvt_range()`](https://wiesehahn.github.io/rvtr/reference/rvt_range.md).

- invert:

  flip the layer after stretching, so high reads dark

## Value

an `rvt_stack`

## See also

[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)
to add layers,
[`rvt_render()`](https://wiesehahn.github.io/rvtr/reference/rvt_render.md)
to write the result

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_hillshade() |> rvt_stack(range = c(0, 1))
#> <rvt_stack> 1 layer, bottom to top
#>   1. file224224db48c0.tif         (background) opacity 1.00  range 0-1
#> Renders at 1000 x 1000, 1 m. Not a raster yet - call rvt_render() to write it.
```
