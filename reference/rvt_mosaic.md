# Present one or many rasters as a single dataset

Passes a single path straight through. For several paths, builds a GDAL
virtual raster (VRT) over them and returns its path - a small XML file,
no pixel data is copied. All functions in this package accept either
form, so calling this directly is only necessary if you want to keep and
reuse the VRT.

## Usage

``` r
rvt_mosaic(dem, vrt_path = fs::file_temp(ext = "vrt"), cl_arg = NULL)
```

## Arguments

- dem:

  a path, or a character vector of paths to mosaic. A directory, or a
  vector of length 1 ending in `.vrt`, is returned unchanged.

- vrt_path:

  where to write the VRT (default: a temp file)

- cl_arg:

  extra arguments passed to `gdalbuildvrt`, e.g.
  `c("-resolution", "highest")`

## Value

path to a single dataset usable by
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md),
[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md),
[`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md).
Pipe-friendly: `tiles |> rvt_mosaic() |> rvt_vat()`.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
rvt_mosaic(dem)  # single path passes through unchanged
#> /home/runner/work/_temp/Library/rvtr/extdata/dtm1.tif
```
