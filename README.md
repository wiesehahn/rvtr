# rvtr

Sky-view factor, positive openness and the archaeological **VAT** blend for R —
following the Relief Visualization Toolbox (Kokalj & Somrak 2019), but with no
Python dependency and no WhiteboxTools.

Built for large inputs: a fused, OpenMP-parallel C++ horizon kernel, tile-by-tile
processing so memory stays bounded, and mosaics of many rasters handled through a
GDAL virtual raster so tiles never produce seams.

> The package is named `rvtr` rather than `relief_vizualization_toolbox` because R
> package names cannot contain underscores. The repository directory keeps the
> longer name.

## Install

Open `relief_vizualization_toolbox.Rproj` in RStudio, then in the console:

```r
install.packages(".", repos = NULL, type = "source")
```

Needs a C++ compiler to build (on Windows, install Rtools first — RStudio
prompts for this automatically if it's missing). Everything else, including
GDAL via `gdalraster`, installs as a regular dependency.

## Use

Runs as-is against the 1 m sample tile bundled with the package. Output path is
optional — omit it and you get a temp file; every function returns the path it
wrote:

```r
library(rvtr)

dem <- system.file("extdata", "dtm1.tif", package = "rvtr")

rvt_svf(dem)       # sky-view factor, 0-1
rvt_openness(dem)  # positive openness, degrees
rvt_vat(dem)       # archaeological VAT (combined), 0-1
```

For your own data, swap `dem` for a path to your DEM, and pass an output path
if you want to keep the file, e.g. `rvt_svf(dem, "svf.tif")`. The DEM is always
the first argument, so these also work with the pipe: `dem |> rvt_vat("vat.tif")`.

### Mosaics

Pass a vector of paths. They are wrapped in a VRT, so the horizon search reads
across file boundaries and the result is identical to processing the same area as
one raster (verified bit-for-bit against a single-file run):

```r
tiles <- list.files("dtm_tiles", "\\.tif$", full.names = TRUE)
tiles |> rvt_vat("vat_mosaic.tif")
```

### Large rasters

Processing is tiled automatically; only one tile is in memory at a time, so input
size is limited by disk rather than RAM. Override the tile edge if you want:

```r
rvt_svf("huge_dem.tif", "svf.tif", tile_size = 4096, progress = TRUE)
```

Outputs are written as **Cloud-Optimized GeoTIFFs** — internally tiled, with
overviews built in. That means reading a result back doesn't have to mean
reading all of it: a viewer (including `rvt_plot()`) can ask for a fast
low-resolution overview instead of decoding the full raster, or a full-resolution
window over just part of a huge one — either way, GDAL only touches the bytes
that window actually needs.

### Threads

Defaults to all cores. `rvt_threads(8)` or `options(rvtr.threads = 8)` to change.

### Plotting

`rvt_plot()` gives a quick look at a result without a separate GIS viewer.
Downsampling happens on read via GDAL, so even a huge raster renders fast:

```r
dem |> rvt_vat("vat.tif") |> rvt_plot()
```

## Correctness

This isn't aiming for bit-for-bit parity with rvt-py — where rvt-py has a bug,
this package does the mathematically correct thing instead (see below).
`testthat::test_local()` checks, using the `dtm1.tif` tile bundled in
`inst/extdata/`:

- sky-view factor and positive openness closely track rvt-py's output
- a mosaic of tiles gives bit-identical results to the same area as one raster
- tile size has no effect on output

That gives confidence the algorithms are implemented correctly without locking
behaviour to another implementation's quirks. The rvt-py comparison outputs
(`inst/extdata/rvt_*.tif`) are **not** committed to git — they're a local
development aid for spot-checking against rvt-py, not a shipped fixture, so
those specific checks skip themselves when the files aren't present. See
`inst/extdata/README.md` for how to regenerate them.

Two rvt-py quirks are deliberately **not** reproduced by default:

- rvt-py treats "no horizon found in this direction" (only reachable beside
  NoData) as a slope of -1000 rather than -∞, off by a fraction of a degree in
  positive openness near NoData holes. This package uses the exact value.
- rvt-py's `blend_overlay()` mutates its `background` argument in place, so the
  VAT openness layer's declared 50% opacity ends up having no effect on its
  output. This package applies the opacity as configured. Pass
  `rvt_vat(rvt_compat = TRUE)` to reproduce rvt-py's output instead, e.g. to
  check results against it directly.

RVT's other padding choice is kept because it's the right one regardless of RVT:
reflect-pad for the horizon search (extrapolates terrain more plausibly at the
raster edge), edge-replicate for the slope/aspect derivatives (standard for
finite differences).
