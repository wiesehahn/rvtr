# rvtr: Relief Visualization for Large Rasters and Mosaics

Topographic visualizations following the Relief Visualization Toolbox
(Kokalj & Somrak 2019): hillshade, multi-directional hillshade, cast
shadow, slope, isotropic and anisotropic sky-view factor, sky
illumination, positive and negative openness, local dominance, simple
and multi-scale relief models, multi-scale topographic position, and the
archaeological Visualization for Archaeological Topography (VAT) blend;
together with aspect, curvature, Laplacian-of-Gaussian edge enhancement
and the geomorphons landform classification. Horizon-search and
box-filter kernels run as fused, OpenMP-parallel C++, and rasters are
processed tile by tile so memory use stays bounded regardless of input
size. A mosaic of many rasters is handled by wrapping it in a GDAL
virtual raster, so searches read across file boundaries seamlessly and
tile edges show no seams.

## Details

Every function takes a DEM path (or, for
[`rvt_mosaic()`](https://wiesehahn.github.io/rvtr/reference/rvt_mosaic.md)
and the functions that call it internally, a vector of paths) as its
first argument, so they compose naturally with the native pipe:

    dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
    dem |> rvt_vat("vat.tif")

    tiles <- list.files("dtm_tiles", "\\.tif$", full.names = TRUE)
    tiles |> rvt_mosaic() |> rvt_svf("svf.tif")

    dem |> rvt_vat("vat.tif") |> rvt_plot()

Shaded relief:
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md),
[`rvt_multi_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_multi_hillshade.md),
[`rvt_shadow()`](https://wiesehahn.github.io/rvtr/reference/rvt_shadow.md).

Openness and sky:
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md),
[`rvt_asvf()`](https://wiesehahn.github.io/rvtr/reference/rvt_asvf.md),
[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md),
[`rvt_openness_negative()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness_negative.md),
[`rvt_sky_illumination()`](https://wiesehahn.github.io/rvtr/reference/rvt_sky_illumination.md),
[`rvt_local_dominance()`](https://wiesehahn.github.io/rvtr/reference/rvt_local_dominance.md).

Sun and energy:
[`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md)
(hours of direct sun),
[`rvt_insolation()`](https://wiesehahn.github.io/rvtr/reference/rvt_insolation.md)
(clear-sky energy, kWh/m2),
[`rvt_shadow()`](https://wiesehahn.github.io/rvtr/reference/rvt_shadow.md).

Surface form:
[`rvt_slope()`](https://wiesehahn.github.io/rvtr/reference/rvt_slope.md),
[`rvt_aspect()`](https://wiesehahn.github.io/rvtr/reference/rvt_aspect.md),
[`rvt_curvature()`](https://wiesehahn.github.io/rvtr/reference/rvt_curvature.md),
[`rvt_log()`](https://wiesehahn.github.io/rvtr/reference/rvt_log.md),
[`rvt_geomorphons()`](https://wiesehahn.github.io/rvtr/reference/rvt_geomorphons.md).

Trend removal:
[`rvt_slrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)
(also
[`rvt_tpi()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)),
[`rvt_dev()`](https://wiesehahn.github.io/rvtr/reference/rvt_dev.md),
[`rvt_msrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_msrm.md),
[`rvt_mstp()`](https://wiesehahn.github.io/rvtr/reference/rvt_mstp.md).

Composite:
[`rvt_vat()`](https://wiesehahn.github.io/rvtr/reference/rvt_vat.md)
ready-made;
[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)
to compose your own;
[`rvt_relight()`](https://wiesehahn.github.io/rvtr/reference/rvt_relight.md)
to light an orthophoto with its terrain, in named
[rvt_relight_styles](https://wiesehahn.github.io/rvtr/reference/rvt_relight_styles.md).

Preparing a surface:
[`rvt_fill()`](https://wiesehahn.github.io/rvtr/reference/rvt_fill.md)
closes NoData gaps,
[`rvt_smooth()`](https://wiesehahn.github.io/rvtr/reference/rvt_smooth.md)
removes noise without rounding off breaks of slope. Both are worth
considering before the noise-sensitive metrics - curvature, LoG,
geomorphons - and in that order, since smoothing cannot average across a
hole:

    dem |> rvt_fill() |> rvt_smooth() |> rvt_curvature()

Support:
[`rvt_mosaic()`](https://wiesehahn.github.io/rvtr/reference/rvt_mosaic.md),
[`rvt_resample()`](https://wiesehahn.github.io/rvtr/reference/rvt_resample.md),
[`rvt_threads()`](https://wiesehahn.github.io/rvtr/reference/rvt_threads.md),
[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md),
[`rvt_palettes()`](https://wiesehahn.github.io/rvtr/reference/rvt_palettes.md)
to browse the palettes
[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)
takes by name.

Data:
[`rvt_data_mapterhorn()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_mapterhorn.md)
streams terrain for anywhere in the world, at the best resolution
published there
([`rvt_data_mapterhorn_sources()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_mapterhorn_sources.md)
lists the sources).
[`rvt_data_lgln()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln.md)
streams terrain, surface and orthophoto data for any place in Lower
Saxony and any flight year
([`rvt_data_lgln_years()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln_years.md)
lists them).

Loading the package sets GDAL's `GTIFF_SRS_SOURCE` option to `"EPSG"`,
unless you have already set it. Many European survey GeoTIFFs describe
their coordinate system in a way that differs trivially from the EPSG
registry, and without this GDAL warns about it every time such a file is
opened. The results are unaffected. To keep GDAL's default behaviour,
set the option to `"GEOKEYS"` yourself, before or after loading.

Not reimplemented here, because `gdalraster` already does them well and
at comparable speed: terrain ruggedness index and roughness, via
`gdalraster::dem_proc(mode = "TRI")` and `mode = "roughness"`. Its
`mode = "TPI"` is fixed to a 3x3 window;
[`rvt_tpi()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)
takes any radius.

## Combining visualizations

No single visualization shows everything, so the usual answer is to lay
several over one another.
[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)
does that with the same blend modes QGIS offers, reading bottom to top:

    dem |> rvt_hillshade() |>
      rvt_blend(rvt_openness(dem), "overlay",  opacity = 0.5) |>
      rvt_blend(rvt_svf(dem),      "multiply", opacity = 0.25) |>
      rvt_plot()

Metrics are computed once and blended as often as you like, so tuning an
opacity costs seconds rather than recomputing a horizon search. For the
soft, evenly-lit look of a rendered terrain model, start from sky
illumination rather than a single sun:

    dem |> rvt_sky_illumination() |>
      rvt_blend(rvt_svf(dem),       "multiply", opacity = 0.3) |>
      rvt_blend(rvt_hillshade(dem), "overlay",  opacity = 0.4) |>
      rvt_plot()

## Choosing a visualization

Kokalj and Hesse's *Guide to Good Practice* (2017) is the standard
reference on this, and its advice is worth summarising here. The
headline point is that **no single visualization is enough**: in almost
every case several are needed, and the choice depends on the topography,
on the size and shape of what you are looking for, and on whether you
are trying to *spot* something or to *interpret* something already
spotted.

Start with
[`rvt_hillshade()`](https://wiesehahn.github.io/rvtr/reference/rvt_hillshade.md)
for an overview, since it gives the most natural-looking read of the
terrain and helps you judge what else is worth running. Then work
outwards roughly in this order, adjusting for terrain:

- flat ground and gentle slopes:

  hillshade with a low sun (under 10 degrees), then
  [`rvt_slrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)
  (~20 m filter),
  [`rvt_local_dominance()`](https://wiesehahn.github.io/rvtr/reference/rvt_local_dominance.md)
  (10-20 m), then
  [`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
  (10 m). Sky-view factor is weakest here - on flat ground it mostly
  shows negative features and is very noise-sensitive.

- moderate slopes:

  hillshade (~30 degrees),
  [`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md)
  (10 m),
  [`rvt_slrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md)
  (~20 m),
  [`rvt_local_dominance()`](https://wiesehahn.github.io/rvtr/reference/rvt_local_dominance.md)
  (10-20 m),
  [`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
  (10 m).

- steep or complex topography:

  hillshade (45 degrees or more),
  [`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md)
  (10 m),
  [`rvt_local_dominance()`](https://wiesehahn.github.io/rvtr/reference/rvt_local_dominance.md)
  (10-20 m),
  [`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md)
  (10 m). Openness is particularly at home here, since removing the
  overall landform means slopes don't saturate the image.

Those radii are in **metres**, and so is every distance this package
takes, so they can be used as published. Each function's own help page
carries the detailed settings and display stretches for its metric.

## Working at more than one scale

The same ground answers differently at different resolutions:
metre-scale detail on a 0.25 m DTM, landform structure on a 2 m version
of it. Since detail cannot be invented, that means coarsening -
[`rvt_resample()`](https://wiesehahn.github.io/rvtr/reference/rvt_resample.md)
does it, and pipes straight into any metric:

    dem |> rvt_resample(2) |> rvt_curvature() |> rvt_plot()

For metrics that take a radius, widening the radius is often the better
way to change scale - it keeps the fine detail in the input instead of
discarding it, and for
[`rvt_slrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_slrm.md),
[`rvt_msrm()`](https://wiesehahn.github.io/rvtr/reference/rvt_msrm.md),
[`rvt_dev()`](https://wiesehahn.github.io/rvtr/reference/rvt_dev.md) and
[`rvt_mstp()`](https://wiesehahn.github.io/rvtr/reference/rvt_mstp.md)
it costs nothing extra.
[`rvt_curvature()`](https://wiesehahn.github.io/rvtr/reference/rvt_curvature.md)
and [`rvt_log()`](https://wiesehahn.github.io/rvtr/reference/rvt_log.md)
have a scale too, but pay for it. Coarsening is the only scale control
where there is no radius at all
([`rvt_slope()`](https://wiesehahn.github.io/rvtr/reference/rvt_slope.md),
[`rvt_aspect()`](https://wiesehahn.github.io/rvtr/reference/rvt_aspect.md)).
See
[`rvt_resample()`](https://wiesehahn.github.io/rvtr/reference/rvt_resample.md)
for the full comparison.

The horizon-search metrics get both at once: `reach` (see
[rvt_reach](https://wiesehahn.github.io/rvtr/reference/rvt_reach.md))
keeps near terrain at full resolution while reading distant terrain from
coarser copies, so microrelief and mountains are accounted for in the
same pass and the cost grows only with the logarithm of the distance.

## References

Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
Visualization: A Guide to Good Practice*. Prostor, kraj, čas 14.
Ljubljana: Založba ZRC.
[doi:10.3986/9789612549848](https://doi.org/10.3986/9789612549848)

## See also

Useful links:

- <https://github.com/wiesehahn/rvtr>

- <https://wiesehahn.github.io/rvtr/>

- Report bugs at <https://github.com/wiesehahn/rvtr/issues>

## Author

**Maintainer**: Jens Wiesehahn <wiesehahn.jens@gmail.com>
