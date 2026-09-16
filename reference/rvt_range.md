# Stretch range for a layer

The `c(lo, hi)` a layer is stretched over before blending, read from the
raster itself. Evaluate it when you build the stack, so the numbers are
fixed before any tile is touched.

## Usage

``` r
rvt_range(path, pct = NULL, band = 1L, max_dim = 1000L)
```

## Arguments

- path:

  raster to measure

- pct:

  optional `c(low, high)` percentages to cut from each tail, e.g.
  `c(2, 98)`. `NULL` (default) uses the full minimum and maximum.

- band:

  band to measure (default 1)

- max_dim:

  when `pct` is given, the longest side of the decimated read the
  percentiles are taken from (default 1000)

## Value

`c(lo, hi)`

## Details

That is the whole point: a range derived per tile would give the same
ground different values depending on which tile it landed in, which is
the one thing this package will not do. Passing a frozen pair sidesteps
the question entirely.

## See also

[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
svf <- rvt_svf(dem)
rvt_range(svf)
#> [1] 0.5444252 1.0020018
rvt_range(svf, pct = c(2, 98))
#> [1] 0.8334779 0.9884582
```
