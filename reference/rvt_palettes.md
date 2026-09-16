# Browse the colour palettes

Draws every palette
[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)
accepts by name as a colour ramp, grouped by kind, so you can pick one
by eye rather than by guessing a name. Pass the name you like straight
to
[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md):
`rvt_plot(svf, "Viridis")`.

## Usage

``` r
rvt_palettes(type = NULL, n = 64)
```

## Arguments

- type:

  `NULL` (default) for all of them, or one of `"sequential"`,
  `"diverging"`, `"divergingx"`, `"qualitative"`

- n:

  number of colours drawn in each ramp (default 64)

## Value

the palette names shown, invisibly

## Details

Which kind to reach for depends on what the values mean:

- **sequential** for magnitudes running one way - sky-view factor,
  slope, openness, elevation. Most relief visualizations want one of
  these, and a plain `"Grays"` is hard to beat for reading shape.

- **diverging** for signed results centred on zero - curvature,
  [`rvt_log()`](https://wiesehahn.github.io/rvtr/reference/rvt_log.md),
  [`rvt_slrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md),
  [`rvt_msrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_msrm.md),
  [`rvt_dev()`](https://wiesehahn.github.io/rvtr/reference/rvt_dev.md) -
  so convex and concave read as opposites. Keep the stretch symmetric
  about zero, as
  [`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)
  explains, or zero will not land on the neutral middle.

- **divergingx** join two sequential ramps with no requirement that they
  balance. Many are familiar from other tools (`"RdBu"`, `"PuOr"`,
  `"Spectral"`), but one arm can end much darker than the other, so
  equal positive and negative values need not look equally strong:
  `"PuOr"` ends at lightness 10 on its purple side against 31 on its
  orange side. `"Cividis"` is filed here too, although it runs dark to
  light and is really sequential. Prefer **diverging** when that balance
  matters.

- **qualitative** for classes with no order, such as
  [`rvt_geomorphons()`](https://wiesehahn.github.io/rvtr/reference/rvt_geomorphons.md):
  a ramp would imply an order the classes do not have.

## See also

[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)

## Examples

``` r
rvt_palettes("diverging")

names <- rvt_palettes("sequential")
```
