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
## would still be ~49 GB. bDOM is heavier again - uncompressed Float32 at 0.2 m,
## 146 MB for one 1 km tile.

## `product` names index this list, and .lgln_product() below detects "the
## caller passed nothing" by comparing against names(), so the formal default
## of rvt_data_lgln() must list these in exactly this order.
.lgln_products <- list(
  dtm = list(host = "dgm", collection = "dgm1", asset = "dgm1-tif", res = 1,
             type = "Float32", label = "DGM1 digital terrain model"),
  dsm = list(host = "dom", collection = "dom1", asset = "dom1-tif", res = 1,
             type = "Float32", label = "DOM1 digital surface model"),
  bdom = list(host = "bdom", collection = "BDOM", asset = "bdom20", res = 0.2,
              type = "Float32", label = "bDOM20 image-based surface model"),
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
    stop('`product` must be one of "dtm", "dsm", "bdom" or "rgb".', call. = FALSE)
  product
}

## Download estimate in MB per km² - the measured figures above. Terrain never
## drops below ~1 MB/km², because there is no overview coarser than 2 m.
.lgln_mb_per_km2 <- function(product, res) {
  if (product == "bdom") {
    # uncompressed Float32, so the cost is exactly the pixels GDAL reads times
    # four bytes, at whichever level it picks: native plus four overviews
    lev <- c(0.2, 0.4, 0.8, 1.6, 3.2)
    mb <- c(100, 25, 6.25, 1.56, 0.39)
    return(mb[max(1L, sum(lev <= res))])
  }
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

.lgln_tile_of <- function(xy) {
  x0 <- floor(xy[1] / 1000) * 1000; y0 <- floor(xy[2] / 1000) * 1000
  c(x0, y0, x0 + 1000, y0 + 1000)
}

## Any location -> list(extent in EPSG:25832, res, kind).
##
## `res` is what the *caller* asked for and may be NULL; `default` is the
## product's own resolution, used only when nothing was asked and the location
## does not bring a grid of its own. The two must stay apart: passing the
## product default as `res` made a raster `x` of any other resolution an error,
## even though the caller had asked for nothing.
.lgln_area <- function(x, res, default = NULL) {
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

  # every other location has no grid of its own, so it is snapped to `res`
  if (is.null(res)) res <- default

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
      return(list(extent = .data_snap(.lgln_tile_of(xy), res), res = res,
                  kind = "point"))
    }
    return(list(extent = .data_snap(as.numeric(sf::st_bbox(g)), res), res = res,
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
      return(list(extent = .data_snap(.lgln_tile_of(xy[1, ]), res), res = res,
                  kind = "point"))
    ext <- c(min(xy[, 1]), min(xy[, 2]), max(xy[, 1]), max(xy[, 2]))
    return(list(extent = .data_snap(ext, res), res = res, kind = "extent"))
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

## What fraction of the extent a set of tiles covers, exactly.
##
## The tiles are axis-aligned rectangles on a whole-kilometre grid, so the area
## of their union is a coordinate sweep: cut the extent at every tile edge, and
## every cell of the resulting grid is either wholly inside a tile or wholly
## outside one. Summing the inside ones is exact, and stays exact if two tiles
## overlap, which adding their areas would not. This replaced a probe grid
## sampled every 500 m, which could not see a narrower gap and could not give
## a fraction.
.lgln_coverage <- function(items, extent) {
  area <- (extent[3] - extent[1]) * (extent[4] - extent[2])
  if (is.null(items) || !nrow(items) || area <= 0) return(0)
  # one row per distinct tile: a year can hold the same tile under several dates
  it <- items[!duplicated(paste(items$xmin, items$ymin, items$xmax, items$ymax)),
              , drop = FALSE]
  brk <- function(v, lo, hi) sort(unique(pmin(pmax(c(lo, hi, v), lo), hi)))
  xs <- brk(c(it$xmin, it$xmax), extent[1], extent[3])
  ys <- brk(c(it$ymin, it$ymax), extent[2], extent[4])
  if (length(xs) < 2L || length(ys) < 2L) return(0)
  covered <- 0
  for (i in seq_len(length(xs) - 1L)) {
    cx <- (xs[i] + xs[i + 1L]) / 2
    inx <- cx > it$xmin & cx < it$xmax
    if (!any(inx)) next
    for (j in seq_len(length(ys) - 1L)) {
      cy <- (ys[j] + ys[j + 1L]) / 2
      if (any(inx & cy > it$ymin & cy < it$ymax))
        covered <- covered + (xs[i + 1L] - xs[i]) * (ys[j + 1L] - ys[j])
    }
  }
  min(1, covered / area)
}

.lgln_covers <- function(items, extent) {
  .lgln_coverage(items, extent) >= 1 - 1e-9
}

## Which items to use, and how much of the extent they cover.
##
## Two different questions, and they were conflated once: *which* survey a tile
## comes from, and whether the result has a hole in it.
##
## `year = NULL` asks for the best available picture of the ground, so every
## tile takes its own latest flight, whatever year that falls in - a partial
## 2024 over a full 2022 gives 2024 where it was flown and 2022 for the rest.
## That is the only way to get both the newest data and full coverage, and it
## is why this does not choose a year at all. The cost is that neighbouring
## tiles can be surveys years apart, so .lgln_get() says so when it happens.
##
## Naming a `year` means that year and no other - one survey, a hole where it
## did not reach. `partial` then decides between that hole and an error, and
## for `year = NULL` it decides the same for ground no year has ever covered.
.lgln_pick <- function(items, extent, year, product, partial = TRUE) {
  if (is.null(items) || !nrow(items))
    stop(sprintf(paste("No LGLN %s data covers this location. rvt_data_lgln() only",
                       "covers Lower Saxony."), product), call. = FALSE)
  years <- sort(unique(items$year), decreasing = TRUE)
  cov <- vapply(years, function(y)
    .lgln_coverage(items[items$year == y, , drop = FALSE], extent), numeric(1))
  full <- years[cov >= 1 - 1e-9]

  if (is.null(year)) {
    it <- items                         # every year, newest per tile below
  } else if (!year %in% years) {
    stop(sprintf("No %s from %s at this location. Years present: %s.", product,
                 year, paste(sort(years), collapse = ", ")), call. = FALSE)
  } else {
    it <- items[items$year == year, , drop = FALSE]
  }

  coverage <- .lgln_coverage(it, extent)
  if (!partial && coverage < 1 - 1e-9) {
    if (is.null(year))
      stop(sprintf(paste("No LGLN %s reaches all of this area - every year",
                         "together covers %s of it, so it lies at the edge of",
                         "the data. Years present: %s. Pass partial = TRUE to",
                         "take what there is, NoData elsewhere."), product,
                   .lgln_pct(coverage), paste(sort(years), collapse = ", ")),
           call. = FALSE)
    stop(sprintf(paste("No %s from %s covers this area (%s of it). Years that do:",
                       "%s. Pass partial = TRUE for the part it does cover,",
                       "NoData elsewhere, or leave `year` unset to fill the rest",
                       "from other years."), product, year, .lgln_pct(coverage),
                 if (length(full)) paste(sort(full), collapse = ", ") else "none"),
         call. = FALSE)
  }

  # newest first, then one row per tile: each tile keeps its latest flight
  it <- it[order(it$date, decreasing = TRUE), , drop = FALSE]
  it <- it[!duplicated(paste(it$xmin, it$ymin, it$xmax, it$ymax)), , drop = FALSE]
  used <- sort(unique(it$year), decreasing = TRUE)
  list(year = .lgln_label(used), years = used, items = it, coverage = coverage,
       full = sort(full))
}

## How to name the survey in a filename, in metadata and in messages: the year,
## or the span when tiles come from several.
.lgln_label <- function(years) {
  if (length(years) == 1L) as.character(years)
  else paste0(min(years), "-", max(years))
}

## A coverage fraction as a percentage, never rounded to "100%" while a gap
## remains - "100% of this area, the rest is NoData" would read as nonsense.
.lgln_pct <- function(x) {
  if (x >= 1 - 1e-9) return("100%")
  sprintf("%.4g%%", min(99.9, max(0.1, 100 * x)))
}

.lgln_guard <- function(product, extent, res, max_mb) {
  km2 <- (extent[3] - extent[1]) * (extent[4] - extent[2]) / 1e6
  mb <- km2 * .lgln_mb_per_km2(product, res)
  if (mb <= max_mb) return(invisible(mb))
  why <- if (product == "bdom" && res < 0.4)
    paste(" The image-based surface model is uncompressed 0.2 m Float32 - 100 MB",
          "per km² at native resolution, about 20 times the orthophoto - so pass",
          "`res` unless you really need every pixel.")
  else if (product %in% c("dtm", "dsm") && res > 2)
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
                      max_mb, threads, partial = TRUE) {
  p <- .lgln_products[[product]]
  .lgln_guard(product, extent, res, max_mb)
  picked <- .lgln_pick(.lgln_search(product, extent), extent, year, product,
                       partial)
  # said before the cache check, so a cached result says the same as a fresh one
  if (length(picked$years) > 1L) {
    n <- table(picked$items$year)[as.character(picked$years)]
    message(sprintf(paste("LGLN %s takes each tile's latest flight, so this",
                          "area combines %d years: %s. Neighbouring tiles can",
                          "be surveys years apart; name a `year` for one",
                          "survey throughout."), product, length(picked$years),
                    paste(sprintf("%d (%d %s)", picked$years, n,
                                  ifelse(n == 1L, "tile", "tiles")),
                          collapse = ", ")))
  }
  if (picked$coverage < 1 - 1e-9)
    message(sprintf("LGLN %s %s covers %s of this area; the rest is %s.%s",
                    product, picked$year, .lgln_pct(picked$coverage),
                    if (product == "rgb") "black, as DOP20 declares no NoData"
                    else "NoData",
                    if (!is.null(year) && length(picked$full))
                      sprintf(" Years covering it all: %s.",
                              paste(picked$full, collapse = ", ")) else ""))

  if (download) {
    if (is.null(out_path))
      out_path <- fs::path(tools::R_user_dir("rvtr", "cache"), "lgln",
                           sprintf("%s_%s_%sm_%s.tif", product, picked$year,
                                   format(res, scientific = FALSE),
                                   paste(format(extent, scientific = FALSE),
                                         collapse = "_")))
    out_path <- .out_path(out_path)
    if (!refresh && fs::file_exists(out_path)) return(out_path)
  }

  if (is.null(.lgln_state$credited)) {
    message("LGLN open data: ", .lgln_credit(), ", ", .lgln_licence,
            ". Add \", Daten geändert\" when you publish anything derived from it.")
    .lgln_state$credited <- TRUE
  }

  vrt <- tryCatch(.lgln_vrt(picked$items$href, extent, res, p$res),
                  error = function(e) stop(sprintf(paste("Could not open the LGLN",
                    "%s files for %s:\n%s"), product, picked$year, conditionMessage(e)),
                    call. = FALSE))
  if (!download) return(.as_path(vrt))

  predictor <- if (p$type == "Byte") "STANDARD" else "FLOATING_POINT"
  md <- c(SOURCE = paste0("LGLN ", p$label, ", ", picked$year),
          ATTRIBUTION = .lgln_credit(), LICENSE = .lgln_licence)
  if (!.data_write_cog(vrt, out_path, predictor, md, threads))
    stop(sprintf(paste("Could not download the LGLN %s for %s. The data is fetched",
                       "on demand, so this needs an internet connection."),
                 product, picked$year), call. = FALSE)
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
#' | `"dsm"` | DOM1 surface model, from lidar: ground plus trees and buildings | 1 m, 1 km tiles | usually one |
#' | `"bdom"` | bDOM20 surface model, matched from the aerial photos | 0.2 m, 1 km tiles | 2021 onwards |
#' | `"rgb"` | DOP20 orthophoto | 0.2 m, 2 km tiles | several, e.g. 2013-2025 |
#'
#' Everything comes in ETRS89 / UTM zone 32N (EPSG:25832), the data's own
#' coordinate system; locations are transformed to find the data, but the data
#' is never reprojected. `rvt_data_lgln_years()` lists which years exist for a
#' place without fetching any data.
#'
#' @section Two surface models, and why the flight date matters:
#' `"dsm"` is measured by lidar; `"bdom"` is computed by matching features
#' between overlapping aerial photographs. bDOM is five times finer - 0.2 m
#' against 1 m - and renders roofs, walls, hedges and open ground more crisply,
#' which is what [rvt_relight()] and the surface metrics want in a town.
#'
#' **Matching only works where the photographs show matchable texture, so over
#' woodland the result depends on when the photographs were taken.** A canopy
#' in leaf gives plenty of texture and is captured well. A bare broadleaf
#' canopy gives almost none, and the match then falls through to whatever does
#' match, usually the forest floor - so bDOM can report *no trees where trees
#' stand*. That failure does not announce itself: the wooded ground comes back
#' smooth and plausible, and a hillshade of it shows forest tracks and gullies
#' as though the canopy were not there.
#'
#' So check the season before trusting bDOM under trees.
#' `rvt_data_lgln_years()` gives the flight dates of every year available, and
#' `year` picks one: a summer flight may describe the canopy well, while the
#' early-spring flights that aerial survey usually favours will not. Where
#' canopy matters regardless - canopy height, forestry, relighting a wooded
#' scene - `"dsm"` is measured rather than inferred, and does not depend on the
#' season.
#'
#' bDOM is also the heaviest product here by a distance: uncompressed 0.2 m
#' Float32, about 100 MB per km² at native resolution against the orthophoto's
#' 4.6 MB. Pass `res` unless you need every pixel; the cost guard stops the call
#' before it starts if you forget.
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
#' Any one year is flown over only part of Lower Saxony, so a single year's
#' tiles may fill your area or stop partway across it.
#'
#' **`year = NULL` (the default) gives every tile its own latest flight**, so
#' you get the newest picture of the ground *and* the widest coverage: where
#' 2024 reaches you get 2024, and where it does not you get whatever year last
#' covered that tile. The cost is that neighbouring tiles can be surveys years
#' apart - trees grown, buildings put up - which can show as a seam. A message
#' names the years whenever more than one is used.
#'
#' **Naming a `year` means that year and no other**, one survey throughout.
#' Where it was not flown you get NoData, and a message names the share
#' covered; [rvt_data_lgln_years()] reports the same share up front, in its
#' `coverage` column.
#'
#' `partial = FALSE` refuses a hole rather than returning one: it errors on a
#' named year that does not reach everywhere, and on ground that no year has
#' ever covered. Within a year, each tile still takes its latest flight -
#' one year's flights are spread over several dates between neighbouring
#' tiles, so years are chosen, not dates.
#'
#' Orthophotos are the exception to the NoData rule. DOP20 declares no NoData
#' value, so the uncovered part of an `"rgb"` result comes back **black**
#' rather than flagged, which no downstream function can tell from a very dark
#' pixel. Use `partial = FALSE` there, or mask the result yourself.
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
#' @param product `"dtm"` (default), `"dsm"`, `"bdom"` or `"rgb"`
#' @param year one flight year, giving a single survey throughout, or `NULL`
#'   (default) to give each tile its own latest flight - see the Years
#'   section, and [rvt_data_lgln_years()] for what exists
#' @param res output resolution in metres, or `NULL` (default) for native -
#'   or, when `x` is a raster, for that raster's own resolution
#' @param download `FALSE` (default) returns a virtual raster reading the
#'   remote files on demand; `TRUE` stores a local copy in the user cache
#' @param partial `TRUE` (default) returns data wherever the chosen tiles
#'   reach and leaves the rest NoData; `FALSE` errors instead of returning a
#'   hole
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
#'   rvt_plot(range = c(0, 255))
#' }
#' @export
rvt_data_lgln <- function(x, product = c("dtm", "dsm", "bdom", "rgb"), year = NULL,
                      res = NULL, download = FALSE, partial = TRUE, max_mb = 500,
                      refresh = FALSE, threads = rvt_threads()) {
  product <- .lgln_product(product)
  if (!is.null(res) && (!is.numeric(res) || length(res) != 1L || !is.finite(res) || res <= 0))
    stop("`res` must be a single positive number of metres, or NULL for native.",
         call. = FALSE)
  if (!is.null(year) && (!is.numeric(year) || length(year) != 1L))
    stop("`year` must be a single year, or NULL for the most recent.", call. = FALSE)
  area <- .lgln_area(x, res, .lgln_products[[product]]$res)
  .lgln_get(area$extent, product, year, area$res, isTRUE(download), NULL,
            refresh, max_mb, threads, isTRUE(partial))
}

#' Years of Lower Saxony data available for a place
#'
#' Lists the flight years LGLN holds for a location, without fetching any
#' data: one row per product and year, with the flight dates, the number of
#' tiles involved, and how much of the area asked for that year's tiles fill.
#'
#' @section What `coverage` and `covers` mean:
#' The catalogue is a patchwork of 1 km and 2 km tiles, and any one year is
#' flown over only part of Lower Saxony. `coverage` is the exact fraction of
#' your area that year's tiles fill, and `covers` is the same thing as a
#' yes-or-no (`coverage` of 1). Anything below 1 means an area lying across
#' the boundary between two flight campaigns, where that year's tiles stop
#' partway across it.
#'
#' It is what to read before naming a `year` in [rvt_data_lgln()], because a
#' named year is used on its own: a year at 0.43 gives you that 43% and
#' NoData for the rest (or, under `partial = FALSE`, an error). Leaving `year`
#' unset instead gives each tile its own latest flight, which covers as much
#' ground as all these rows together - so a set of rows that are each partial
#' can still add up to a complete raster, assembled from several years.
#'
#' A point is always answered with the single 1 km tile containing it, so
#' `coverage` is 1 for every year listed. Only an extent, an sf bbox or a
#' raster grid can lie across a seam.
#'
#' @inheritParams rvt_data_lgln
#' @param product one or more of `"dtm"`, `"dsm"`, `"bdom"`, `"rgb"` (default all)
#' @return a data frame, one row per product and year:
#'   \describe{
#'     \item{`product`}{`"dtm"`, `"dsm"`, `"bdom"` or `"rgb"`}
#'     \item{`year`}{the flight year}
#'     \item{`first_date`, `last_date`}{earliest and latest flight date among
#'       that year's tiles, as `"YYYY-MM-DD"`. They differ when neighbouring
#'       tiles were flown on different days of the same campaign}
#'     \item{`tiles`}{how many distinct tiles that year contributes to the
#'       area}
#'     \item{`coverage`}{the fraction of the area those tiles fill, 0 to 1 -
#'       see above. Below 1, naming that year in [rvt_data_lgln()] returns
#'       that share and leaves the rest NoData}
#'     \item{`covers`}{`coverage == 1`, as a convenience: the years that can
#'       be named on their own and still fill the whole area}
#'   }
#' @seealso [rvt_data_lgln()]
#' @examples
#' \donttest{
#' rvt_data_lgln_years(c(9.9464, 51.6317))
#' }
#' @export
rvt_data_lgln_years <- function(x, product = c("dtm", "dsm", "bdom", "rgb")) {
  if (!is.character(product) || !length(product) ||
      !all(product %in% names(.lgln_products)))
    stop('`product` must be any of "dtm", "dsm", "bdom" and "rgb".', call. = FALSE)
  out <- lapply(product, function(pr) {
    area <- .lgln_area(x, NULL, .lgln_products[[pr]]$res)
    items <- .lgln_search(pr, area$extent)
    if (is.null(items) || !nrow(items)) return(NULL)
    do.call(rbind, lapply(sort(unique(items$year)), function(y) {
      it <- items[items$year == y, , drop = FALSE]
      cov <- .lgln_coverage(it, area$extent)
      data.frame(product = pr, year = y, first_date = min(it$date),
                 last_date = max(it$date),
                 tiles = length(unique(paste(it$xmin, it$ymin, it$xmax, it$ymax))),
                 coverage = round(cov, 4), covers = cov >= 1 - 1e-9,
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
