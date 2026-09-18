# Laplacian of Gaussian (edge enhancement)

Finds edges at a chosen scale: smooth the terrain by a Gaussian of width
`sigma`, then take its second derivative. The result is near zero on
smooth ground and swings sharply positive and negative on either side of
a break of slope, so banks, ditch lips, terrace edges and wall lines
come out as crisp paired lines.

## Usage

``` r
rvt_log(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  sigma = 2,
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

- sigma:

  standard deviation of the Gaussian smoothing, in **map units**
  (metres, normally). Unlike the window radii elsewhere it is not
  snapped to whole cells - it is a continuous scale - but it cannot be
  smaller than one cell (default 2)

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

The reason to prefer this over plain
[`rvt_curvature()`](https://wiesehahn.github.io/rvtr/reference/rvt_curvature.md)
is the `sigma`: it sets which size of edge you are looking for and
suppresses everything finer, so a noisy lidar DTM gives a usable image
instead of static.

Kokalj and Hesse suggest using it less as a standalone product than as
an **overlay**: averaged into a sky-view factor or local dominance image
it sharpens edges and, usefully, counteracts the saturation those
metrics suffer on steep slopes when the stretch has been set for gentle
ground. Values are centred on zero and **convex ground is positive**, as
in
[`rvt_curvature()`](https://wiesehahn.github.io/rvtr/reference/rvt_curvature.md)
and the relief models: the lip of a bank or the top of an edge reads
positive, a ditch bottom or the foot of a slope negative. That is the
negative of the raw Laplacian, so where Kokalj and Hesse recommend an
inverted greyscale for the Laplacian, a plain greyscale gives the same
picture here. Display it with a symmetric stretch and a diverging
palette.

## How it works

The Gaussian is separable, so smoothing is a horizontal pass followed by
a vertical one - cost grows with `sigma`, not with `sigma` squared. The
kernel is truncated at three standard deviations, past which the weights
are below a thousandth of the peak. A 3x3 Laplacian then takes the
second derivative of the smoothed surface.

NoData cells are dropped from the smoothing and the weights
renormalised, so a hole blurs across rather than bleeding emptiness
outwards, but the hole itself stays NoData. Overlap is the Gaussian
radius plus one.

## Tuning

- `sigma` is in **map units** and is the whole story: it sets the size
  of edge detected. Around 1-2 picks out sharp, narrow features; 4-8
  finds broader terraces and scarps while ignoring texture. Kokalj and
  Hesse work with a Laplacian radius of about 3 cells.

- Larger `sigma` costs proportionally more and, more importantly, blurs
  away the fine detail - it is a choice about scale, not quality.

## References

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

## See also

[`rvt_curvature()`](https://wiesehahn.github.io/rvtr/reference/rvt_curvature.md)
for the unsmoothed version.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_log(sigma = 3)
```
