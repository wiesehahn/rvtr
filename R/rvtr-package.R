#' @details
#' Every function takes a DEM path (or, for [rvt_mosaic()] and the functions
#' that call it internally, a vector of paths) as its first argument, so they
#' compose naturally with the native pipe:
#'
#' ```r
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#' dem |> rvt_vat("vat.tif")
#'
#' tiles <- list.files("dtm_tiles", "\\.tif$", full.names = TRUE)
#' tiles |> rvt_mosaic() |> rvt_svf("svf.tif")
#'
#' dem |> rvt_vat("vat.tif") |> rvt_plot()
#' ```
#'
#' Shaded relief: [rvt_hillshade()], [rvt_multi_hillshade()], [rvt_shadow()].
#'
#' Openness and sky: [rvt_svf()], [rvt_asvf()], [rvt_openness()],
#' [rvt_openness_negative()], [rvt_sky_illumination()],
#' [rvt_local_dominance()].
#'
#' Sun and energy: [rvt_daylight()] (hours of direct sun),
#' [rvt_insolation()] (clear-sky energy, kWh/m2), [rvt_shadow()].
#'
#' Surface form: [rvt_slope()], [rvt_aspect()], [rvt_curvature()],
#' [rvt_log()], [rvt_geomorphons()].
#'
#' Trend removal: [rvt_slrm()] (also [rvt_tpi()]), [rvt_dev()], [rvt_msrm()],
#' [rvt_mstp()].
#'
#' Composite: [rvt_vat()] ready-made; [rvt_blend()] to compose your own.
#'
#' Preparing a surface: [rvt_fill()] closes NoData gaps, [rvt_smooth()]
#' removes noise without rounding off breaks of slope. Both are worth
#' considering before the noise-sensitive metrics - curvature, LoG,
#' geomorphons - and in that order, since smoothing cannot average across a
#' hole:
#'
#' ```r
#' dem |> rvt_fill() |> rvt_smooth() |> rvt_curvature()
#' ```
#'
#' Support: [rvt_mosaic()], [rvt_resample()], [rvt_threads()], [rvt_plot()],
#' [rvt_palettes()] to browse the palettes `rvt_plot()` takes by name.
#'
#' Data: [rvt_data_lgln()] streams terrain, surface and orthophoto data for any place
#' in Lower Saxony and any flight year ([rvt_data_lgln_years()] lists them).
#'
#' Loading the package sets GDAL's `GTIFF_SRS_SOURCE` option to `"EPSG"`,
#' unless you have already set it. Many European survey GeoTIFFs describe their
#' coordinate system in a way that differs trivially from the EPSG registry,
#' and without this GDAL warns about it every time such a file is opened. The
#' results are unaffected. To keep GDAL's default behaviour, set the option to
#' `"GEOKEYS"` yourself, before or after loading.
#'
#' Not reimplemented here, because `gdalraster` already does them well and at
#' comparable speed: terrain ruggedness index and roughness, via
#' `gdalraster::dem_proc(mode = "TRI")` and `mode = "roughness"`. Its
#' `mode = "TPI"` is fixed to a 3x3 window; [rvt_tpi()] takes any radius.
#'
#' @section Combining visualizations:
#' No single visualization shows everything, so the usual answer is to lay
#' several over one another. [rvt_blend()] does that with the same blend modes
#' QGIS offers, reading bottom to top:
#'
#' ```r
#' dem |> rvt_hillshade() |>
#'   rvt_blend(rvt_openness(dem), "overlay",  opacity = 0.5) |>
#'   rvt_blend(rvt_svf(dem),      "multiply", opacity = 0.25) |>
#'   rvt_plot()
#' ```
#'
#' Metrics are computed once and blended as often as you like, so tuning an
#' opacity costs seconds rather than recomputing a horizon search. For the soft,
#' evenly-lit look of a rendered terrain model, start from sky illumination
#' rather than a single sun:
#'
#' ```r
#' dem |> rvt_sky_illumination() |>
#'   rvt_blend(rvt_svf(dem),       "multiply", opacity = 0.3) |>
#'   rvt_blend(rvt_hillshade(dem), "overlay",  opacity = 0.4) |>
#'   rvt_plot()
#' ```
#'
#' @section Choosing a visualization:
#' Kokalj and Hesse's *Guide to Good Practice* (2017) is the standard
#' reference on this, and its advice is worth summarising here. The headline
#' point is that **no single visualization is enough**: in almost every case
#' several are needed, and the choice depends on the topography, on the size
#' and shape of what you are looking for, and on whether you are trying to
#' *spot* something or to *interpret* something already spotted.
#'
#' Start with [rvt_hillshade()] for an overview, since it gives the most
#' natural-looking read of the terrain and helps you judge what else is worth
#' running. Then work outwards roughly in this order, adjusting for terrain:
#'
#' \describe{
#'   \item{flat ground and gentle slopes}{hillshade with a low sun (under 10
#'     degrees), then [rvt_slrm()] (~20 m filter), [rvt_local_dominance()]
#'     (10-20 m), then [rvt_openness()] (10 m). Sky-view factor is weakest
#'     here - on flat ground it mostly shows negative features and is very
#'     noise-sensitive.}
#'   \item{moderate slopes}{hillshade (~30 degrees), [rvt_svf()] (10 m),
#'     [rvt_slrm()] (~20 m), [rvt_local_dominance()] (10-20 m),
#'     [rvt_openness()] (10 m).}
#'   \item{steep or complex topography}{hillshade (45 degrees or more),
#'     [rvt_svf()] (10 m), [rvt_local_dominance()] (10-20 m), [rvt_openness()]
#'     (10 m). Openness is particularly at home here, since removing the
#'     overall landform means slopes don't saturate the image.}
#' }
#'
#' Those radii are in **metres**, and so is every distance this package takes,
#' so they can be used as published. Each function's own help page carries the
#' detailed settings and display stretches for its metric.
#'
#' @section Working at more than one scale:
#' The same ground answers differently at different resolutions: metre-scale
#' detail on a 0.25 m DTM, landform structure on a 2 m version of it. Since
#' detail cannot be invented, that means coarsening - [rvt_resample()] does it,
#' and pipes straight into any metric:
#'
#' ```r
#' dem |> rvt_resample(2) |> rvt_curvature() |> rvt_plot()
#' ```
#'
#' For metrics that take a radius, widening the radius is often the better way
#' to change scale - it keeps the fine detail in the input instead of
#' discarding it, and for [rvt_slrm()], [rvt_msrm()], [rvt_dev()] and
#' [rvt_mstp()] it costs nothing extra. [rvt_curvature()] and [rvt_log()] have
#' a scale too, but pay for it. Coarsening is the only scale control where
#' there is no radius at all ([rvt_slope()], [rvt_aspect()]). See
#' [rvt_resample()] for the full comparison.
#'
#' The horizon-search metrics get both at once: `reach` (see [rvt_reach])
#' keeps near terrain at full resolution while reading distant terrain from
#' coarser copies, so microrelief and mountains are accounted for in the same
#' pass and the cost grows only with the logarithm of the distance.
#'
#' @references
#' Kokalj, Ž. and Hesse, R. (2017) *Airborne Laser Scanning Raster Data
#' Visualization: A Guide to Good Practice*. Prostor, kraj, čas 14. Ljubljana:
#' Založba ZRC. \doi{10.3986/9789612549848}
#'
#' @keywords internal
#' @useDynLib rvtr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
