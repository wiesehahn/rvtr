# Fill gaps in a raster

Interpolates across NoData holes so the surface is continuous. Worth
doing before anything else, because a hole is not simply "no answer
here": it spreads. The derivative metrics
([`rvt_slope()`](https://wiesehahn.github.io/rvtr/reference/rvt_slope.md),
[`rvt_curvature()`](https://wiesehahn.github.io/rvtr/reference/rvt_curvature.md))
lose a ring of cells around every hole, and the horizon metrics treat a
hole as *transparent* - light passes straight through ground that was
never measured - so a DSM with gaps in its canopy will read as brighter
than the real thing.

## Usage

``` r
rvt_fill(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  max_distance = 100,
  smooth_iterations = 0,
  threads = rvt_threads(),
  overwrite = FALSE
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

- max_distance:

  how far to search for valid cells, in **map units** (metres, normally;
  default 100)

- smooth_iterations:

  3x3 smoothing passes applied to the filled cells only (default 0)

- threads:

  C++ threads (see
  [`rvt_threads()`](https://wiesehahn.github.io/rvtr/reference/rvt_threads.md))

- overwrite:

  recompute even if `out_path` exists

## Value

`out_path`, invisibly

## How it works

GDAL's `GDALFillNodata()` (via
[`gdalraster::fillNodata()`](https://firelab.github.io/gdalraster/reference/fillNodata.html)):
for each cell in a hole, the value is interpolated from the nearest
valid cells in each of four directions, weighted by inverse distance,
followed by optional smoothing passes over the interpolated values only.

Cells that already had data are never touched - verified bit for bit -
so filling can only add information, never alter what was measured.

Unlike the metrics here this is **not tiled**: GDAL fills the band as a
whole. On a raster far larger than memory, fill the tiles individually
before mosaicking.

## Tuning

- `max_distance` is how far the interpolation will reach to find valid
  ground, in **map units**. Holes wider than twice this are left partly
  unfilled, which is usually what you want: a 2 km gap in coverage
  should not be invented from its edges.

- `smooth_iterations` runs 3x3 averaging passes over the filled cells,
  softening the spokes that inverse-distance interpolation can leave in
  a large hole. 0 is fine for scattered single-cell gaps; 2-3 helps on
  anything larger.

## See also

[`rvt_smooth()`](https://wiesehahn.github.io/rvtr/reference/rvt_smooth.md),
which removes noise rather than gaps - run this first, since smoothing
cannot sensibly cross a hole.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_fill()
```
