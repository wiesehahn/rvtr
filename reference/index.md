# Package index

## Overview

- [`rvtr`](https://wiesehahn.github.io/rvtr/reference/rvtr-package.md)
  [`rvtr-package`](https://wiesehahn.github.io/rvtr/reference/rvtr-package.md)
  : rvtr: Relief Visualization for Large Rasters and Mosaics
- [`rvt_reach`](https://wiesehahn.github.io/rvtr/reference/rvt_reach.md)
  : Reaching further with a multi-resolution search

## Shading and sun

- [`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)
  : Hillshade
- [`rvt_multi_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_multi_hillshade.md)
  : Multi-directional hillshade
- [`rvt_shadow()`](https://wiesehahn.github.io/rvtr/reference/rvt_shadow.md)
  : Cast shadow
- [`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md)
  : Time in daylight
- [`rvt_insolation()`](https://wiesehahn.github.io/rvtr/reference/rvt_insolation.md)
  : Clear-sky solar irradiation

## Sky and horizon

- [`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md) :
  Sky-view factor
- [`rvt_asvf()`](https://wiesehahn.github.io/rvtr/reference/rvt_asvf.md)
  : Anisotropic sky-view factor
- [`rvt_sky_illumination()`](https://wiesehahn.github.io/rvtr/reference/rvt_sky_illumination.md)
  : Sky illumination
- [`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
  : Positive openness
- [`rvt_openness_negative()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness_negative.md)
  : Negative openness

## Local relief

- [`rvt_slrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)
  [`rvt_tpi()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)
  : Simple local relief model (SLRM)
- [`rvt_msrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_msrm.md)
  : Multi-scale relief model (MSRM)
- [`rvt_dev()`](https://wiesehahn.github.io/rvtr/reference/rvt_dev.md) :
  Deviation from mean elevation (DEV)
- [`rvt_local_dominance()`](https://wiesehahn.github.io/rvtr/reference/rvt_local_dominance.md)
  : Local dominance
- [`rvt_mstp()`](https://wiesehahn.github.io/rvtr/reference/rvt_mstp.md)
  : Multi-scale topographic position (MSTP)

## Surface shape

- [`rvt_slope()`](https://wiesehahn.github.io/rvtr/reference/rvt_slope.md)
  : Slope
- [`rvt_aspect()`](https://wiesehahn.github.io/rvtr/reference/rvt_aspect.md)
  : Aspect
- [`rvt_curvature()`](https://wiesehahn.github.io/rvtr/reference/rvt_curvature.md)
  : Curvature
- [`rvt_log()`](https://wiesehahn.github.io/rvtr/reference/rvt_log.md) :
  Laplacian of Gaussian (edge enhancement)
- [`rvt_geomorphons()`](https://wiesehahn.github.io/rvtr/reference/rvt_geomorphons.md)
  : Geomorphons
- [`rvt_geomorphon_classes`](https://wiesehahn.github.io/rvtr/reference/rvt_geomorphon_classes.md)
  : Geomorphon landform classes

## Preparation

- [`rvt_fill()`](https://wiesehahn.github.io/rvtr/reference/rvt_fill.md)
  : Fill gaps in a raster
- [`rvt_smooth()`](https://wiesehahn.github.io/rvtr/reference/rvt_smooth.md)
  : Feature-preserving smoothing
- [`rvt_resample()`](https://wiesehahn.github.io/rvtr/reference/rvt_resample.md)
  : Resample a raster to a coarser resolution
- [`rvt_mosaic()`](https://wiesehahn.github.io/rvtr/reference/rvt_mosaic.md)
  : Present one or many rasters as a single dataset

## Composites and blending

- [`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md) :
  Archaeological VAT (combined)
- [`rvt_preset_general`](https://wiesehahn.github.io/rvtr/reference/rvt_presets.md)
  [`rvt_preset_flat`](https://wiesehahn.github.io/rvtr/reference/rvt_presets.md)
  : VAT terrain presets
- [`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)
  : Blend a layer over a raster
- [`rvt_blend_modes`](https://wiesehahn.github.io/rvtr/reference/rvt_blend_modes.md)
  : Blend modes
- [`rvt_stack()`](https://wiesehahn.github.io/rvtr/reference/rvt_stack.md)
  : Start a blend stack
- [`rvt_range()`](https://wiesehahn.github.io/rvtr/reference/rvt_range.md)
  : Stretch range for a layer
- [`rvt_render()`](https://wiesehahn.github.io/rvtr/reference/rvt_render.md)
  : Render a blend stack to a raster
- [`rvt_relight()`](https://wiesehahn.github.io/rvtr/reference/rvt_relight.md)
  : Relight an orthophoto with its terrain
- [`rvt_relight_styles`](https://wiesehahn.github.io/rvtr/reference/rvt_relight_styles.md)
  : Lighting styles for rvt_relight()

## Plotting

- [`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)
  : Quick-look plot of a raster
- [`rvt_plot3d()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot3d.md)
  : Interactive 3D view of a terrain model
- [`rvt_palettes()`](https://wiesehahn.github.io/rvtr/reference/rvt_palettes.md)
  : Browse the colour palettes

## Data

- [`rvt_data_mapterhorn()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_mapterhorn.md)
  : Terrain for anywhere in the world, from Mapterhorn
- [`rvt_data_mapterhorn_sources()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_mapterhorn_sources.md)
  : Which Mapterhorn sources cover a place
- [`rvt_data_lgln()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln.md)
  : Lower Saxony elevation and orthophotos, for any place
- [`rvt_data_lgln_years()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln_years.md)
  : Years of Lower Saxony data available for a place

## Settings

- [`rvt_threads()`](https://wiesehahn.github.io/rvtr/reference/rvt_threads.md)
  : Number of threads used by the C++ kernels
