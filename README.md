# rvtr

Relief visualization and terrain analysis for R. Sky-view factor, openness,
local dominance, relief models, solar radiation, curvature, geomorphons and
more, from a DEM straight to a Cloud-Optimized GeoTIFF.

- **Large rasters**: processed tile by tile, so memory use stays bounded and
  input size is limited by disk rather than RAM.
- **Seamless mosaics**: pass a vector of files and the result is identical to
  processing the same area as one raster, with no seams at file boundaries.
- **Fast**: fused, multi-threaded C++ kernels.
- **Pipeable**: every function takes a raster path first and returns the path
  it wrote.

The visualizations follow the
[Relief Visualization Toolbox](https://github.com/EarthObservation/RVT_py)
(Kokalj & Somrak 2019), reimplemented independently in R and C++.

## Installation

```r
# install.packages("pak")
pak::pak("wiesehahn/rvtr")
```

Building from source needs a C++ compiler (on Windows, install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)). GDAL comes with the
`gdalraster` dependency. Install `rgl` too if you want 3D views.

## Quick start

```r
library(rvtr)

dem <- system.file("extdata", "dtm1.tif", package = "rvtr")  # 1 m sample tile

dem |> rvt_svf() |> rvt_plot()        # temporary file, plotted
rvt_svf(dem, "svf.tif")               # or keep the result
dem |> rvt_vat() |> rvt_plot()        # ready-made composite for spotting earthworks
```

## What's included

| | functions | |
|---|---|---|
| **Shading and sun** | `rvt_hillshade()`, `rvt_multi_hillshade()`, `rvt_shadow()` | shaded relief from one or several directions; true cast shadow |
| | `rvt_daylight()`, `rvt_insolation()` | hours of direct sun; clear-sky solar energy in kWh/m² |
| **Sky and horizon** | `rvt_svf()`, `rvt_asvf()`, `rvt_sky_illumination()` | how much sky each cell sees, so nothing hides in shadow |
| | `rvt_openness()`, `rvt_openness_negative()` | banks and mounds, or ditches and hollow ways, brought forward |
| **Local relief** | `rvt_slrm()` / `rvt_tpi()`, `rvt_msrm()`, `rvt_dev()` | height above or below the smoothed terrain, at one scale or many |
| | `rvt_local_dominance()`, `rvt_mstp()` | subtle, spread-out features; RGB composite by scale |
| **Surface shape** | `rvt_slope()`, `rvt_aspect()`, `rvt_curvature()`, `rvt_log()` | steepness, orientation, 14 curvature types, edges at a chosen scale |
| | `rvt_geomorphons()` | a landform class per cell: ridge, spur, hollow, valley… |
| **Preparation** | `rvt_fill()`, `rvt_smooth()`, `rvt_resample()`, `rvt_mosaic()` | fill NoData, feature-preserving smoothing, coarsen, combine files |
| **Display** | `rvt_vat()`, `rvt_blend()`, `rvt_render()` | composites and 13 blend modes |
| | `rvt_relight()` | orthophotos lit by their own terrain, with cast and soft shadows |
| | `rvt_plot()`, `rvt_plot3d()`, `rvt_palettes()` | quick 2D and interactive 3D views |
| **Data** | `rvt_data_lgln()` | open DTM, DSM and orthophotos for Lower Saxony |

Every help page explains how the method works, how to tune it, and which
settings suit flat, moderate and steep terrain. `?rvtr` gives an overview of
which visualization to use where.

## Good to know

**Distances are in map units**, normally metres, never pixels. `reach = 10`
looks 10 m out whether the DEM has 1 m or 0.25 m cells, so values from the
literature can be used as they are.

**Mosaics**: pass several files at once.

```r
tiles <- list.files("dtm_tiles", "\\.tif$", full.names = TRUE)
tiles |> rvt_svf("svf.tif", progress = TRUE)
```

**Long reach stays cheap.** The horizon metrics read distant terrain from
coarser copies of the DEM, so the cost grows with the logarithm of the distance
and accuracy is hardly affected. See `?rvt_reach`.

```r
dem |> rvt_openness(reach = 400)
```

**Preparation chains**: fill holes, then smooth, then analyse, optionally at a
coarser scale.

```r
dem |> rvt_fill() |> rvt_smooth() |> rvt_resample(2) |> rvt_curvature() |> rvt_plot()
```

**Threads**: all cores by default; change with `rvt_threads(8)`.

## Blending

Layer finished rasters with the blend modes a GIS offers. The piped raster is
the background, and nothing is computed until you plot or render:

```r
dem |> rvt_hillshade() |>
  rvt_blend(rvt_openness(dem), "overlay",  opacity = 0.5) |>
  rvt_blend(rvt_svf(dem),      "multiply", opacity = 0.25) |>
  rvt_render("relief.tif")
```

The [blending vignette](vignettes/blending.Rmd) walks through all thirteen
modes.

## Sample data

`rvt_data_lgln()` streams open data from Lower Saxony's survey authority (LGLN)
for any point or extent in the state:

```r
pt <- c(9.9464, 51.6317)   # longitude, latitude
rvt_data_lgln(pt, "dtm") |> rvt_hillshade() |> rvt_plot()
```

Data licence CC BY 4.0: © GeoBasis-DE/LGLN (year of download), plus
"Daten geändert" for derived products.

## References

- Kokalj, Ž. & Somrak, M. (2019). Why not a single image? Combining
  visualizations to facilitate fieldwork and on-screen mapping. *Remote
  Sensing* 11(7), 747. <https://doi.org/10.3390/rs11070747>
- Kokalj, Ž. & Hesse, R. (2017). *Airborne Laser Scanning Raster Data
  Visualization: A Guide to Good Practice*. Založba ZRC.
  <https://doi.org/10.3986/9789612549848>

Individual help pages cite the sources for each method.

## License

MIT
