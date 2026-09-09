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
#' Main functions: [rvt_svf()], [rvt_openness()], [rvt_vat()], [rvt_mosaic()],
#' [rvt_threads()], [rvt_plot()].
#'
#' @keywords internal
#' @useDynLib rvtr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
