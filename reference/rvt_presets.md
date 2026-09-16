# VAT terrain presets

The two parameter sets that
[`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md)
renders and averages. `rvt_preset_general` suits ordinary relief;
`rvt_preset_flat` looks further out and filters noise harder, to bring
up the subtle micro-relief typical of flat ground.

## Usage

``` r
rvt_preset_general

rvt_preset_flat
```

## Format

A list with:

- `sun_elevation`:

  height of the sun above the horizon, in degrees, for the hillshade
  layer. Lower angles throw longer shadows and exaggerate faint relief -
  hence 15 for flat terrain against 35 for general.

- `slope`:

  the slope range, in degrees, stretched across the layer's full
  brightness range. Slopes past the upper end are clipped, so a narrower
  range gives more contrast among gentle slopes (15 for flat terrain, 50
  for general).

- `svf`:

  the same idea for the sky-view factor layer - the sky-fraction range
  that gets stretched. `c(0.9, 1)` on flat terrain spreads out what
  would otherwise be an almost invisible band of values.

- `opns`:

  the same again for the positive openness layer, in degrees.

- `reach`:

  how far the horizon search behind the openness and sky-view factor
  layers looks, in map units (metres, normally).

- `noise_removal`:

  0-3, ignoring progressively more of the innermost part of each search
  ray.

An object of class `list` of length 6.

An object of class `list` of length 6.

## Details

Both are plain lists, so the way to adjust one is a modified copy:
`modifyList(rvt_preset_general, list(reach = 20))`.

## See also

[`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md)
