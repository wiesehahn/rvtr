# Deviation from mean elevation (DEV)

How far a cell sits above or below its surroundings, measured in local
standard deviations rather than metres. Positive on local highs,
negative in local lows, and roughly 0 wherever the ground follows its
own trend.

## Usage

``` r
rvt_dev(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  radius = 20,
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
  mosaics are read across file boundaries, so tiles do not produce seams

- out_path:

  output GeoTIFF path (default: a temp file)

- radius:

  radius of the window a cell is compared against, in **map units**
  (metres, normally) (default 20)

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

This is
[`rvt_slrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)
divided by the local roughness, and that division is the whole point: it
makes the number comparable between smooth and rugged ground. A 20 cm
bank in a flat ploughed field and a 2 m terrace on a broken hillside can
both come out around 2, because each is measured against how much its
own neighbourhood varies. Values are dimensionless, and behave much like
a z-score - most ground falls within about -2 to 2.

The trade is that it says nothing about actual height: a strong DEV on
very flat ground may be a few centimetres of noise. Use
[`rvt_slrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)
when you want metres and DEV when you want prominence.

## How it works

For a square window of `radius`, `(cell - mean) / sd`, with mean and
standard deviation both taken over the window from summed-area tables -
so the cost is independent of the radius. The variance is computed about
the window mean rather than by differencing sums of squares, which
matters on flat ground at high elevation (see
[`rvt_mstp()`](https://wiesehahn.github.io/rvtr/reference/rvt_mstp.md),
which maximises this same quantity over a range of scales).

Terrain is mirrored at the raster edge rather than smeared outwards, so
a window hanging over the edge still sees varied ground and the standard
deviation doesn't collapse.
[`rvt_slrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)
replicates the edge instead, since it has no divisor to protect - so
within `radius` of the border the two can disagree slightly, including
in sign. Further in than that they agree exactly.

## Tuning

- `radius` is in **map units** and sets the neighbourhood a cell is
  judged against - the size of feature it responds to. Too small and
  everything looks locally average; too large and small features are
  drowned by the surrounding landform.

- If you don't know the right scale,
  [`rvt_mstp()`](https://wiesehahn.github.io/rvtr/reference/rvt_mstp.md)
  tries a whole range and keeps whichever is most pronounced.

## See also

[`rvt_slrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)
for the same thing in elevation units,
[`rvt_mstp()`](https://wiesehahn.github.io/rvtr/reference/rvt_mstp.md)
for the multi-scale version.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
rvt_dev(dem)
```
