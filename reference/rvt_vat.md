# Archaeological VAT (combined)

Visualization for Archaeological Topography: a ready-made composite for
spotting earthworks, combining four different views of the terrain into
a single greyscale image so you don't have to read several rasters side
by side. Values are 0-1, ready to display directly.

## Usage

``` r
rvt_vat(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  num_directions = 16,
  presets = list(rvt_preset_general, rvt_preset_flat),
  tile_size = NULL,
  threads = rvt_threads(),
  rvt_compat = FALSE,
  quality = 90,
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

  output GeoTIFF path (default: a temp file). A folder that does not
  exist yet is created.

- num_directions:

  number of compass directions scanned (default 16)

- presets:

  list of two parameter sets, general first then flat; see
  [rvt_presets](https://wiesehahn.github.io/rvtr/reference/rvt_presets.md)

- tile_size:

  tile edge in pixels; NULL picks a size targeting roughly 256 MB per
  working matrix

- threads:

  C++ threads (see
  [`rvt_threads()`](https://wiesehahn.github.io/rvtr/reference/rvt_threads.md))

- rvt_compat:

  reproduce a bug in rvt-py's blend_overlay() that makes the openness
  layer's declared 50% opacity have no effect. `FALSE` (default) applies
  the opacity as configured; `TRUE` is only useful for comparing output
  against rvt-py.

- quality:

  for a `.webp` or `.jpg` `out_path`, compression quality 1-100 (default
  90)

- overwrite:

  recompute even if `out_path` exists

- progress:

  report per-tile progress

## Value

`out_path`, invisibly. `dem` is the first argument, so this is
pipe-friendly: `dem |> rvt_vat("vat.tif")`.

## Details

Each ingredient contributes something the others lack: hillshade gives
the eye the familiar sense of relief, slope sharpens edges, positive
openness lifts convex features such as banks and mounds, and sky-view
factor darkens enclosed ground such as ditches. The whole thing is
rendered twice - once with settings for ordinary relief, once tuned for
flat ground where features are subtle - and the two are averaged, so a
single run works reasonably across both without retuning.

If you want to interpret actual numbers rather than look at a picture,
this is the wrong tool: it is a display product, and the blending
discards the physical units of its ingredients. Use
[`rvt_msrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_msrm.md)
(heights in metres) or
[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
(angles in degrees) for that.

## How it works

For each of the two presets
([rvt_presets](https://wiesehahn.github.io/rvtr/reference/rvt_presets.md)),
the DEM is turned into four layers, each stretched to 0-1 over a range
set by the preset, then stacked from the bottom up:

- **hillshade** at the preset's sun elevation, from the north-west - the
  background everything else is laid over;

- **slope**, inverted so flat ground is bright, blended over it at half
  strength;

- **positive openness**, blended in Overlay mode, which lightens what is
  already light and darkens what is already dark, so convex features
  gain contrast without washing the image out;

- **sky-view factor**, multiplied in at a quarter strength, which only
  darkens - deepening enclosed ground.

The two resulting images are then averaged. The horizon searches use
each preset's own `reach` and `noise_removal`, so the flat-terrain pass
looks further out and filters more aggressively than the general one.
Computation is at native resolution, tile by tile.

## Tuning

- The two presets are ordinary lists, so the easiest adjustment is a
  modified copy - for example
  `modifyList(rvt_preset_general, list(reach = 20))` to make the general
  pass respond to larger features. See
  [rvt_presets](https://wiesehahn.github.io/rvtr/reference/rvt_presets.md)
  for the fields and what they mean.

- To render a single terrain type rather than the average of two, pass
  the same preset twice, e.g.
  `presets = list(rvt_preset_flat, rvt_preset_flat)`. This is worth
  doing when your area really is uniformly flat, since the
  general-terrain pass otherwise dilutes the contrast the flat preset
  was tuned to give.

- `num_directions` is shared by both passes and behaves as in
  [`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md).

- `rvt_compat` exists only for checking output against rvt-py and should
  be left alone otherwise; see its parameter description.

## Output

A Cloud-Optimized GeoTIFF, values 0-1. Give `out_path` a `.webp` or
`.jpg` extension to write a plain greyscale picture instead, with no
GeoTIFF made along the way and no georeferencing; see
[`rvt_image()`](https://wiesehahn.github.io/rvtr/reference/rvt_image.md)
for the formats.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_vat()
```
