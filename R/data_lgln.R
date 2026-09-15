## Lower Saxony open data from LGLN's STAC APIs, for any place and flight year.
##
## Streamed, not downloaded, by default: the result is a small VRT over the
## remote Cloud-Optimized GeoTIFFs, so GDAL fetches only the blocks a read
## needs, from overviews when coarser. Data is never reprojected - everything
## is native EPSG:25832.
##
## Measured per 1 km tile (GDAL network statistics): orthophoto 4.61 MB at
## 0.2 m but 0.41 MB at 1 m, from overviews; terrain 3.40 MB at 1 m and 1.03 MB
## at 2 m *or any coarser resolution*, because the DGM1/DOM1 COGs have a single
## 2 m overview. That floor is why the cost guard exists: Lower Saxony at 100 m
## would still be ~49 GB.

.lgln_products <- list(
  dtm = list(host = "dgm", collection = "dgm1", asset = "dgm1-tif", res = 1,
             type = "Float32", label = "DGM1 digital terrain model"),
  dsm = list(host = "dom", collection = "dom1", asset = "dom1-tif", res = 1,
             type = "Float32", label = "DOM1 digital surface model"),
  rgb = list(host = "dop", collection = "DOP", asset = "dop20_rgb", res = 0.2,
             type = "Byte", label = "DOP20 RGB orthophoto")
)

.lgln_licence <- "CC BY 4.0"
.lgln_credit <- function(year = format(Sys.Date(), "%Y")) {
  paste0("© GeoBasis-DE/LGLN ", year)
}
.lgln_state <- new.env(parent = emptyenv())

.lgln_product <- function(product) {
  if (identical(product, names(.lgln_products))) return("dtm")    # the default
  # our own message: match.arg()'s arrives translated into the session language
  if (!is.character(product) || length(product) != 1L ||
      !product %in% names(.lgln_products))
    stop('`product` must be one of "dtm", "dsm" or "rgb".', call. = FALSE)
  product
}

## Download estimate in MB per km² - the measured figures above. Terrain never
## drops below ~1 MB/km², because there is no overview coarser than 2 m.
.lgln_mb_per_km2 <- function(product, res) {
  if (product == "rgb") {
    if (res <= 0.2) 4.61
    else if (res < 1) max(0.41, 4.61 * (0.2 / res)^2)
    else max(0.005, 0.41 / res^2)
  } else {
    if (res <= 1) 3.40 else 1.03
  }
}

## Tile extent in EPSG:25832 from an LGLN item id: "..._32_<east km>_<north
## km>_<size km>_ni_..." e.g. dop20rgbi_32_564_5720_2_ni_2025-03-08.
.lgln_tile_extent <- function(id) {
  m <- regmatches(id, regexec("_32_(\\d+)_(\\d+)_(\\d+)_ni_", id))[[1]]
  if (length(m) != 4L) return(rep(NA_real_, 4))
  v <- as.numeric(m[2:4]) * 1000
  c(v[1], v[2], v[1] + v[3], v[2] + v[3])
}

.lgln_snap <- function(extent, res) {
  r <- function(v) round(v, 6)
  c(r(floor(extent[1] / res + 1e-9) * res), r(floor(extent[2] / res + 1e-9) * res),
    r(ceiling(extent[3] / res - 1e-9) * res), r(ceiling(extent[4] / res - 1e-9) * res))
}

.lgln_tile_of <- function(xy) {
  x0 <- floor(xy[1] / 1000) * 1000; y0 <- floor(xy[2] / 1000) * 1000
  c(x0, y0, x0 + 1000, y0 + 1000)
}

## Any location -> list(extent in EPSG:25832, res, kind).
.lgln_area <- function(x, res) {
  # a raster: exactly its grid, which must already be EPSG:25832
  if (is.character(x) && length(x) == 1L && fs::file_exists(x)) {
    ds <- methods::new(gdalraster::GDALRaster, x, read_only = TRUE)
    on.exit(ds$close())
    epsg <- gdalraster::srs_find_epsg(ds$getProjection())
    if (!identical(epsg, "EPSG:25832"))
      stop("A raster given as `x` must be in EPSG:25832, the data's own CRS.",
           call. = FALSE)
    rr <- abs(ds$res()[1])
    if (!is.null(res) && !isTRUE(all.equal(res, rr)))
      stop("`res` can't differ from the resolution of the raster given as `x`.",
           call. = FALSE)
    return(list(extent = as.numeric(ds$bbox()), res = rr, kind = "raster"))
  }

  # sf objects carry their own CRS
  if (inherits(x, c("sf", "sfc", "sfg", "bbox"))) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Install the sf package to pass sf objects as `x`.", call. = FALSE)
    g <- if (inherits(x, "bbox")) sf::st_as_sfc(x) else sf::st_geometry(x)
    if (is.na(sf::st_crs(g)))
      stop("The sf object given as `x` has no CRS.", call. = FALSE)
    g <- sf::st_transform(g, 25832)
    if (length(g) == 1L && all(sf::st_geometry_type(g) == "POINT")) {
      xy <- as.numeric(sf::st_coordinates(g))[1:2]
      return(list(extent = .lgln_snap(.lgln_tile_of(xy), res), res = res,
                  kind = "point"))
    }
    return(list(extent = .lgln_snap(as.numeric(sf::st_bbox(g)), res), res = res,
                kind = "extent"))
  }

  # plain numbers are WGS84 longitude and latitude
  if (is.numeric(x) && length(x) %in% c(2L, 4L)) {
    lon <- if (length(x) == 2L) x[1] else x[c(1, 3, 1, 3)]
    lat <- if (length(x) == 2L) x[2] else x[c(2, 2, 4, 4)]
    if (any(abs(lon) > 180) || any(abs(lat) > 90))
      stop(paste("Plain numbers in `x` are WGS84 longitude and latitude. For",
                 "EPSG:25832 coordinates, pass an sf point or bbox with its CRS."),
           call. = FALSE)
    xy <- gdalraster::transform_xy(cbind(lon, lat), "EPSG:4326", "EPSG:25832")
    if (length(x) == 2L)
      return(list(extent = .lgln_snap(.lgln_tile_of(xy[1, ]), res), res = res,
                  kind = "point"))
    ext <- c(min(xy[, 1]), min(xy[, 2]), max(xy[, 1]), max(xy[, 2]))
    return(list(extent = .lgln_snap(ext, res), res = res, kind = "extent"))
  }

  stop(paste("`x` must be a longitude/latitude point c(lon, lat), a bbox",
             "c(xmin, ymin, xmax, ymax) in longitude/latitude, an sf point or",
             "geometry, or the path to an EPSG:25832 raster."), call. = FALSE)
}

.lgln_lonlat_bbox <- function(extent) {
  ll <- gdalraster::transform_xy(
    cbind(extent[c(1, 3, 1, 3)], extent[c(2, 2, 4, 4)]), "EPSG:25832", "EPSG:4326")
  c(min(ll[, 1]), min(ll[, 2]), max(ll[, 1]), max(ll[, 2]))
}

## Catalogue search -> data frame of usable items intersecting the extent.
.lgln_search <- function(product, extent) {
  p <- .lgln_products[[product]]
  ll <- .lgln_lonlat_bbox(extent)
  url <- sprintf("https://%s.stac.lgln.niedersachsen.de/search?collections=%s&bbox=%s&limit=1000",
                 p$host, p$collection, paste(sprintf("%.6f", ll), collapse = ","))
  feats <- list()
  for (page in seq_len(50L)) {
    res <- tryCatch(jsonlite::fromJSON(url, simplifyVector = FALSE),
                    error = function(e) e)
    if (inherits(res, "error"))
      stop(sprintf("Could not reach the LGLN catalogue at\n  %s\n%s", url,
                   conditionMessage(res)), call. = FALSE)
    feats <- c(feats, res$features)
    nxt <- Filter(function(l) identical(l$rel, "next") &&
                    (is.null(l$method) || identical(toupper(l$method), "GET")),
                  if (is.null(res$links)) list() else res$links)
    if (!length(nxt)) break
    url <- nxt[[1]]$href
  }

  rows <- lapply(feats, function(f) {
    href <- f$assets[[p$asset]]$href
    if (is.null(href)) return(NULL)
    # 2025 is catalogued at 10 cm too, but those files are not published (every
    # link 404s); only 20 cm orthophotos are usable
    if (product == "rgb") {
      px <- f$properties$bodenpixelgroesse
      if (!startsWith(f$id, "dop20") || (!is.null(px) && as.numeric(px) != 20))
        return(NULL)
    }
    te <- .lgln_tile_extent(f$id)
    if (anyNA(te)) return(NULL)
    dt <- substr(f$properties$datetime, 1, 10)
    data.frame(id = f$id, date = dt, year = as.integer(substr(dt, 1, 4)),
               xmin = te[1], ymin = te[2], xmax = te[3], ymax = te[4],
               href = href, stringsAsFactors = FALSE)
  })
  items <- do.call(rbind, rows)
  if (is.null(items)) return(items)
  hit <- items$xmin < extent[3] & items$xmax > extent[1] &
         items$ymin < extent[4] & items$ymax > extent[2]
  items[hit, , drop = FALSE]
}

## Points sampled over the extent, to test whether a year's tiles cover it.
.lgln_probe <- function(extent) {
  step <- max(500, (extent[3] - extent[1]) / 200, (extent[4] - extent[2]) / 200)
  xs <- unique(pmin(pmax(seq(extent[1], extent[3], by = step), extent[1] + 0.01),
                    extent[3] - 0.01))
  ys <- unique(pmin(pmax(seq(extent[2], extent[4], by = step), extent[2] + 0.01),
                    extent[4] - 0.01))
  as.matrix(expand.grid(c(xs, extent[3] - 0.01), c(ys, extent[4] - 0.01)))
}

.lgln_covers <- function(items, extent) {
  pts <- .lgln_probe(extent)
  all(vapply(seq_len(nrow(pts)), function(i)
    any(pts[i, 1] >= items$xmin & pts[i, 1] < items$xmax &
        pts[i, 2] >= items$ymin & pts[i, 2] < items$ymax), logical(1)))
}

## The year to use, and one item (the latest flight in that year) per tile.
.lgln_pick <- function(items, extent, year, product) {
  if (is.null(items) || !nrow(items))
    stop(sprintf(paste("No LGLN %s data covers this location. rvt_data_lgln() only",
                       "covers Lower Saxony."), product), call. = FALSE)
  years <- sort(unique(items$year), decreasing = TRUE)
  full <- years[vapply(years, function(y)
    .lgln_covers(items[items$year == y, , drop = FALSE], extent), logical(1))]
  if (is.null(year)) {
    if (!length(full))
      stop(sprintf(paste("No single year of %s covers this whole area, it lies at",
                         "the edge of the data. Years present: %s."), product,
                   paste(sort(years), collapse = ", ")), call. = FALSE)
    year <- full[1]
  } else if (!year %in% full) {
    stop(sprintf("No %s from %s covers this area. Years that do: %s.", product,
                 year, if (length(full)) paste(sort(full), collapse = ", ") else "none"),
         call. = FALSE)
  }
  it <- items[items$year == year, , drop = FALSE]
  it <- it[order(it$date, decreasing = TRUE), , drop = FALSE]
  it <- it[!duplicated(paste(it$xmin, it$ymin, it$xmax, it$ymax)), , drop = FALSE]
  list(year = as.integer(year), items = it)
}

.lgln_guard <- function(product, extent, res, max_mb) {
  km2 <- (extent[3] - extent[1]) * (extent[4] - extent[2]) / 1e6
  mb <- km2 * .lgln_mb_per_km2(product, res)
  if (mb <= max_mb) return(invisible(mb))
  why <- if (product != "rgb" && res > 2)
    paste(" The terrain and surface files have no overview coarser than 2 m, so",
          "a coarse `res` still reads about 1 MB per km²; use a coarse elevation",
          "source for large areas.")
  else ""
  stop(sprintf(paste("This request would download an estimated %.0f MB (%.0f km²",
                     "of %s at %s m), above `max_mb = %s`.%s Reduce the area, or",
                     "raise `max_mb` if you mean it."),
               mb, km2, product, format(res), format(max_mb), why), call. = FALSE)
}

## Build the lazy raster: a native-resolution crop, or a warped VRT when
## resampling. Both reference the remote COGs, so nothing is downloaded yet.
.lgln_vrt <- function(hrefs, extent, res, native) {
  src <- paste0("/vsicurl/", hrefs)
  te <- sprintf("%.6f", extent)
  tr <- format(res, scientific = FALSE)
  out <- fs::file_temp(ext = "vrt")
  if (isTRUE(all.equal(res, native))) {
    gdalraster::buildVRT(out, src, cl_arg = c("-te", te, "-tr", tr, tr), quiet = TRUE)
  } else {
    mosaic <- fs::file_temp(ext = "vrt")
    gdalraster::buildVRT(mosaic, src, quiet = TRUE)
    gdalraster::warp(mosaic, out, t_srs = "EPSG:25832",
                     cl_arg = c("-of", "VRT", "-te", te, "-tr", tr, tr,
                                "-r", "cubic"),
                     quiet = TRUE)
  }
  out
}

## Fetch one product for an extent: lazy VRT, or a cached COG with download.
.lgln_get <- function(extent, product, year, res, download, out_path, refresh,
                      max_mb, threads) {
  p <- .lgln_products[[product]]
  .lgln_guard(product, extent, res, max_mb)
  picked <- .lgln_pick(.lgln_search(product, extent), extent, year, product)

  if (download) {
    if (is.null(out_path))
      out_path <- fs::path(tools::R_user_dir("rvtr", "cache"), "lgln",
                           sprintf("%s_%d_%sm_%s.tif", product, picked$year,
                                   format(res, scientific = FALSE),
                                   paste(format(extent, scientific = FALSE),
                                         collapse = "_")))
    out_path <- .as_path(out_path)
    if (!refresh && fs::file_exists(out_path)) return(out_path)
  }

  if (is.null(.lgln_state$credited)) {
    message("LGLN open data: ", .lgln_credit(), ", ", .lgln_licence,
            ". Add \", Daten geändert\" when you publish anything derived from it.")
    .lgln_state$credited <- TRUE
  }

  vrt <- tryCatch(.lgln_vrt(picked$items$href, extent, res, p$res),
                  error = function(e) stop(sprintf(paste("Could not open the LGLN",
                    "%s files for %d:\n%s"), product, picked$year, conditionMessage(e)),
                    call. = FALSE))
  if (!download) return(.as_path(vrt))

  fs::dir_create(fs::path_dir(out_path))
  # written to a hidden file beside the target, then moved into place, so an
  # interrupted download never leaves a truncated file that looks valid
  tmp <- fs::path(fs::path_dir(out_path), paste0(".", fs::path_file(out_path), ".part"))
  .rm_path(tmp)
  on.exit(.rm_path(tmp), add = TRUE)
  predictor <- if (p$type == "Byte") "STANDARD" else "FLOATING_POINT"
  ok <- tryCatch({
    gdalraster::translate(vrt, tmp,
      cl_arg = c("-of", "COG",
                 "-co", "COMPRESS=DEFLATE",
                 "-co", paste0("PREDICTOR=", predictor),
                 "-co", "RESAMPLING=CUBIC",
                 "-co", "OVERVIEWS=IGNORE_EXISTING",
                 "-co", paste0("NUM_THREADS=", threads),
                 "-mo", paste0("SOURCE=LGLN ", p$label, ", ", picked$year),
                 "-mo", paste0("ATTRIBUTION=", .lgln_credit()),
                 "-mo", paste0("LICENSE=", .lgln_licence)),
      quiet = TRUE)
    fs::file_exists(tmp)
  }, error = function(e) FALSE)
  if (!isTRUE(ok))
    stop(sprintf(paste("Could not download the LGLN %s for %d. The data is fetched",
                       "on demand, so this needs an internet connection."),
                 product, picked$year), call. = FALSE)
  .rm_path(out_path)
  fs::file_move(tmp, out_path)
  out_path
}

#' Lower Saxony elevation and orthophotos, for any place
#'
#' Terrain model, surface model or orthophoto for any location in Lower Saxony,
#' from the state survey authority (LGLN), in any flight year it holds. Give a
#' point to get the 1 km tile it lies in, an extent to get exactly that area, or
#' a raster to get data on its grid.
#'
#' Nothing is downloaded by default. The result is a small virtual raster that
#' points at LGLN's published Cloud-Optimized GeoTIFFs, and data moves only when
#' something reads it - and only the parts that read needs. Pass it to any
#' function in this package. `download = TRUE` stores a local copy instead, for
#' offline work or repeated sessions.
#'
#' | product | what | native grid | years |
#' |---|---|---|---|
#' | `"dtm"` | DGM1 terrain model, the bare ground | 1 m, 1 km tiles | usually one |
#' | `"dsm"` | DOM1 surface model, ground plus trees and buildings | 1 m, 1 km tiles | usually one |
#' | `"rgb"` | DOP20 orthophoto | 0.2 m, 2 km tiles | several, e.g. 2013-2025 |
#'
#' Everything comes in ETRS89 / UTM zone 32N (EPSG:25832), the data's own
#' coordinate system; locations are transformed to find the data, but the data
#' is never reprojected. `rvt_data_lgln_years()` lists which years exist for a
#' place without fetching any data.
#'
#' @section Location:
#' * **A point** returns the 1 x 1 km tile containing it. Plain numbers
#'   `c(lon, lat)` are WGS84; an sf point uses its own coordinate system, so
#'   EPSG:25832 coordinates work as an sf point.
#' * **An extent** returns exactly that area, snapped outward to whole cells:
#'   an `sf::st_bbox()`, any sf geometry, or plain WGS84
#'   `c(xmin, ymin, xmax, ymax)`.
#' * **A raster path** returns data on exactly that raster's grid, which must
#'   be EPSG:25832.
#'
#' @section Resolution and download size:
#' `res = NULL` keeps native resolution. A number resamples, cubic, onto a grid
#' of whole multiples of `res`, so products fetched for the same place at the
#' same `res` line up exactly and can go straight into [rvt_blend()].
#'
#' Coarser requests read from the files' overviews, which cuts orthophoto
#' downloads roughly tenfold at 1 m. **The terrain and surface files have no
#' overview coarser than 2 m**, though, so even a 100 m request still reads
#' about 1 MB per km² - all of Lower Saxony would be around 50 GB. Before
#' anything is read, the download is estimated and the call stops above
#' `max_mb`. For large areas at coarse resolution, use a coarse elevation
#' source instead.
#'
#' @section Years:
#' `year = NULL` uses the most recent year that covers the whole area. Flights
#' in one year can be spread over several dates between neighbouring tiles, so
#' years are chosen, not dates; within a year each tile uses its latest flight.
#'
#' @section Licence:
#' LGLN open geodata, licensed CC BY 4.0 under section 7 of LGLN's terms of
#' use. Credit them as **© GeoBasis-DE/LGLN** followed by the year you
#' downloaded the data, and add **", Daten geändert"** (data modified) for
#' anything derived from them.
#'
#' @param x location: a WGS84 `c(lon, lat)` point or `c(xmin, ymin, xmax,
#'   ymax)` extent, an sf point, geometry or bbox, or the path to an
#'   EPSG:25832 raster
#' @param product `"dtm"` (default), `"dsm"` or `"rgb"`
#' @param year flight year, or `NULL` (default) for the most recent that covers
#'   the whole area
#' @param res output resolution in metres, or `NULL` (default) for native
#' @param download `FALSE` (default) returns a virtual raster reading the
#'   remote files on demand; `TRUE` stores a local copy in the user cache
#' @param max_mb stop before reading if the estimated download exceeds this
#'   many megabytes (default 500)
#' @param refresh with `download = TRUE`, fetch again even if cached
#' @inheritParams rvt_svf
#' @return path to the raster (a VRT, or a Cloud-Optimized GeoTIFF with
#'   `download = TRUE`), so it pipes into any function
#' @seealso [rvt_data_lgln_years()]
#' @examples
#' \donttest{
#' pt <- c(9.9464, 51.6317)                      # Burg Hardenberg
#' rvt_data_lgln_years(pt, "rgb")
#' rvt_data_lgln(pt, "dtm") |> rvt_hillshade() |> rvt_plot()
#' rvt_data_lgln(pt, "rgb", year = 2013, res = 1) |>
#'   rvt_plot(minmax_def = c(0, 0, 0, 255, 255, 255))
#' }
#' @export
rvt_data_lgln <- function(x, product = c("dtm", "dsm", "rgb"), year = NULL,
                      res = NULL, download = FALSE, max_mb = 500,
                      refresh = FALSE, threads = rvt_threads()) {
  product <- .lgln_product(product)
  if (!is.null(res) && (!is.numeric(res) || length(res) != 1L || !is.finite(res) || res <= 0))
    stop("`res` must be a single positive number of metres, or NULL for native.",
         call. = FALSE)
  if (!is.null(year) && (!is.numeric(year) || length(year) != 1L))
    stop("`year` must be a single year, or NULL for the most recent.", call. = FALSE)
  area <- .lgln_area(x, if (is.null(res)) .lgln_products[[product]]$res else res)
  .lgln_get(area$extent, product, year, area$res, isTRUE(download), NULL,
            refresh, max_mb, threads)
}

#' Years of Lower Saxony data available for a place
#'
#' Lists the flight years LGLN holds for a location, without fetching any
#' data: one row per product and year, with the flight dates, the number of
#' tiles involved, and whether that year covers the whole area.
#'
#' @inheritParams rvt_data_lgln
#' @param product one or more of `"dtm"`, `"dsm"`, `"rgb"` (default all)
#' @return a data frame with columns `product`, `year`, `first_date`,
#'   `last_date`, `tiles` and `covers`
#' @seealso [rvt_data_lgln()]
#' @examples
#' \donttest{
#' rvt_data_lgln_years(c(9.9464, 51.6317))
#' }
#' @export
rvt_data_lgln_years <- function(x, product = c("dtm", "dsm", "rgb")) {
  if (!is.character(product) || !length(product) ||
      !all(product %in% names(.lgln_products)))
    stop('`product` must be any of "dtm", "dsm" and "rgb".', call. = FALSE)
  out <- lapply(product, function(pr) {
    area <- .lgln_area(x, .lgln_products[[pr]]$res)
    items <- .lgln_search(pr, area$extent)
    if (is.null(items) || !nrow(items)) return(NULL)
    do.call(rbind, lapply(sort(unique(items$year)), function(y) {
      it <- items[items$year == y, , drop = FALSE]
      data.frame(product = pr, year = y, first_date = min(it$date),
                 last_date = max(it$date),
                 tiles = length(unique(paste(it$xmin, it$ymin, it$xmax, it$ymax))),
                 covers = .lgln_covers(it, area$extent),
                 stringsAsFactors = FALSE)
    }))
  })
  res <- do.call(rbind, out)
  if (is.null(res))
    stop("No LGLN data covers this location. rvt_data_lgln() only covers Lower Saxony.",
         call. = FALSE)
  rownames(res) <- NULL
  res
}
