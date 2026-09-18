# Multi-scale relief model (MSRM)

Strips away both the broad shape of the landscape and fine-grained
noise, leaving only relief within a chosen size range.

## Usage

``` r
rvt_msrm(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  feature_min = 0,
  feature_max = 20,
  scaling_factor = 2,
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

- feature_min, feature_max:

  smallest/largest feature size to keep, in the DEM's horizontal map
  units, typically metres; `feature_min` is clamped up to the pixel
  resolution if smaller (defaults 0 and 20)

- scaling_factor:

  positive integer spacing the sampled radii
  (`radius = n^scaling_factor`); 1 samples every radius, higher values
  sample fewer, more widely spaced ones (default 2)

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

Values are in the DEM's elevation units and read directly as height
above or below the local surface: positive for anything raised (banks,
mounds, walls), negative for anything cut in (ditches, hollow ways,
pits). A value of 0.3 on a metric DEM means the cell sits about 0.3 m
proud of its surroundings at the scales you selected - which makes this
unusually easy to interpret compared with the angle- and ratio-based
metrics here.

Because the regional trend is removed, a low bank on a steep hillside
shows up just as clearly as one on level ground - the situation where
hillshade and sky-view factor struggle most. **Flat ground gives 0, and
so does a uniform slope**: only departures from the local trend survive.

## How it works

The DEM is smoothed repeatedly, by averaging over square windows of
increasing size, and it is the *differences between* those smoothed
versions that form the result:

- The window radii come from `feature_min` and `feature_max`, converted
  from map units to pixels using the raster's own resolution, and spaced
  apart as `n^scaling_factor`. At the defaults on a 1 m DEM that gives
  radii of 0, 1, 4, 9 and 16 pixels.

- Subtracting each smoothed surface from the previous, less-smoothed one
  leaves exactly the detail that step of smoothing removed - the relief
  at that particular scale.

- Those differences are averaged into a single band-pass image, keeping
  what lies between the smallest and largest scale and discarding both
  finer noise and broader landscape form.

The window averages are computed from a summed-area table, so their cost
does not grow with window size - a large `feature_max` is very nearly as
cheap as a small one. NoData cells are excluded from the averages rather
than counted as zero. Computation is at native resolution, tile by tile,
with tiles overlapping by the largest radius.

## Tuning

- `feature_min` is the smallest thing kept. Leaving it at 0 (it is
  clamped up to one cell) keeps the sharpest detail; raise it to a few
  times the cell size on noisy lidar or vegetation-affected DEMs to
  suppress speckle.

- `feature_max` sets how much landscape form is removed - features much
  larger than it fade out. Lower it to flatten strong topography harder
  and concentrate on small features; raise it to keep broader ones.

- `scaling_factor` controls how many scales are sampled between the two.
  1 samples every radius: the most scales, smoothest response, best
  sensitivity to intermediate feature sizes, slightly slower. The
  default 2 uses fewer, more widely spaced scales, which gives more
  contrast. Higher values are punchier still but can step over
  intermediate sizes.

- `feature_min`/`feature_max` are in **map units**, as every distance in
  this package is, so the same settings mean the same ground distances
  regardless of resolution.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_msrm()
```
