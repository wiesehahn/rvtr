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
#' Metrics: [rvt_hillshade()], [rvt_multi_hillshade()], [rvt_shadow()],
#' [rvt_slope()], [rvt_svf()], [rvt_asvf()], [rvt_sky_illumination()],
#' [rvt_openness()], [rvt_openness_negative()], [rvt_local_dominance()],
#' [rvt_slrm()], [rvt_msrm()], [rvt_mstp()], [rvt_vat()].
#'
#' Support: [rvt_mosaic()], [rvt_threads()], [rvt_plot()].
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
#' All of those radii are in **metres**, while this package takes pixels -
#' multiply by four on a 0.25 m DEM. Each function's own help page carries the
#' detailed settings and display stretches for its metric.
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
