## Global terrain from Mapterhorn (mapterhorn.com), for any place on Earth.
##
## Mapterhorn publishes Terrarium-encoded WebP tiles in Web Mercator: 30 m
## Copernicus GLO-30 everywhere up to zoom 12, and national LiDAR models at
## higher zooms where they exist. GDAL's PMTiles driver reads vector tiles
## only, so the archives are not usable; the XYZ tile server is, through the
## WMS driver's TMS mode.
##
## The chain is entirely virtual until something reads it:
##   WMS tiles (RGB, EPSG:3857) at one zoom
##   -> derived VRT decoding Terrarium to Float32, missing tiles as NoData
##   -> mosaic of several zooms, coarse under fine, so the 30 m base shows
##      through wherever a finer source has no tile
##   -> warp to a metric CRS (local UTM by default) at `res`.
##
## Web Mercator cannot be returned as is: its map units are stretched by
## 1/cos(latitude), and every distance in this package is taken as metres.

.mt_tiles_url <- "https://tiles.mapterhorn.com/${z}/${x}/${y}.webp"
.mt_coverage_url <- "https://single-archive-tiles.mapterhorn.com/coverage/%d/%d/%d.mvt"
.mt_attribution_url <- "https://download.mapterhorn.com/attribution.json"
.mt_merc <- 20037508.342789244     # half the Web Mercator world width, metres
.mt_planet_zoom <- 12L             # the global base exists at every tile to here
.mt_max_zoom <- 21L                # 512 * 2^21 still fits a raster dimension
.mt_state <- new.env(parent = emptyenv())

## Measured: tiles at zoom 12-16 are 45-160 kB (flat lowland to the Alps).
## The guard takes the upper end, so it errs towards stopping.
.mt_mb_per_tile <- 0.15

## Ground size of one tile pixel at zoom z and latitude lat, metres.
.mt_px <- function(z, lat) 2 * .mt_merc / (512 * 2^z) * cos(lat * pi / 180)

## The coarsest zoom whose pixels are no larger than `res` on the ground.
## Mapterhorn builds each source up to exactly this zoom: verified for 1 m
## (z16 at 51.6 and 59.9 N), 0.5 m (z17, 46 N), 0.25 m (z18, 48.8 N), 2 m (z15,
## 42.9 S) and 5 m (z14, 52 N), each with a 404 one zoom further.
.mt_zoom_for <- function(res, lat) {
  z <- ceiling(log2(2 * .mt_merc * cos(lat * pi / 180) / (512 * res)) - 1e-9)
  as.integer(pmin(pmax(z, 0L), .mt_max_zoom))
}

## UTM zone of a WGS84 point, as an EPSG code. Standard zones only - no
## Norway/Svalbard exceptions; pass `crs` for those.
.mt_utm <- function(lon, lat) {
  zone <- min(max(floor((lon + 180) / 6) + 1, 1), 60)
  sprintf("EPSG:%d", (if (lat >= 0) 32600 else 32700) + zone)
}

.mt_crs <- function(crs) {
  if (is.null(crs)) return(NULL)
  if (is.numeric(crs) && length(crs) == 1L) crs <- sprintf("EPSG:%d", as.integer(crs))
  if (!is.character(crs) || length(crs) != 1L)
    stop("`crs` must be an EPSG code such as 25832, a string such as \"EPSG:25832\", or WKT.",
         call. = FALSE)
  # quiet: GDAL's own message for an unknown code would print before ours
  gdalraster::push_error_handler("quiet")
  ok <- tryCatch(nzchar(gdalraster::srs_to_wkt(crs)), error = function(e) FALSE,
                 warning = function(w) FALSE)
  gdalraster::pop_error_handler()
  if (!isTRUE(ok))
    stop(sprintf("`crs` = \"%s\" is not a coordinate system GDAL recognises.", crs),
         call. = FALSE)
  if (!gdalraster::srs_is_projected(crs))
    stop(paste("`crs` must be a projected coordinate system in metres. Distances",
               "in this package are map units, so degrees would not work."),
         call. = FALSE)
  crs
}

## Any location -> list(extent in `crs` (unsnapped), crs, ll = WGS84 bbox,
## res = fixed resolution or NULL, grid = TRUE if the extent is a raster's grid).
.mt_area <- function(x, size, crs) {
  crs <- .mt_crs(crs)
  square <- function(lon, lat) {
    crs <- if (is.null(crs)) .mt_utm(lon, lat) else crs
    xy <- gdalraster::transform_xy(cbind(lon, lat), "EPSG:4326", crs)
    list(extent = c(xy[1] - size / 2, xy[2] - size / 2, xy[1] + size / 2, xy[2] + size / 2),
         crs = crs, point = TRUE)
  }
  box <- function(ll) {
    crs <- if (is.null(crs)) .mt_utm(mean(ll[c(1, 3)]), mean(ll[c(2, 4)])) else crs
    list(extent = gdalraster::bbox_transform(ll, "EPSG:4326", crs), crs = crs)
  }

  a <- if (is.character(x) && length(x) == 1L && fs::file_exists(x)) {
    # a raster: exactly its grid and its own CRS
    ds <- methods::new(gdalraster::GDALRaster, x, read_only = TRUE)
    on.exit(ds$close())
    rcrs <- ds$getProjection()
    if (!nzchar(rcrs) || !gdalraster::srs_is_projected(rcrs))
      stop("A raster given as `x` must be in a projected coordinate system in metres.",
           call. = FALSE)
    if (!is.null(crs) && !gdalraster::srs_is_same(crs, rcrs))
      stop("`crs` can't differ from the coordinate system of the raster given as `x`.",
           call. = FALSE)
    list(extent = as.numeric(ds$bbox()), crs = rcrs, res = abs(ds$res()[1]), grid = TRUE)
  } else if (inherits(x, c("sf", "sfc", "sfg", "bbox"))) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Install the sf package to pass sf objects as `x`.", call. = FALSE)
    g <- if (inherits(x, "bbox")) sf::st_as_sfc(x) else sf::st_geometry(x)
    if (inherits(g, "sfg")) g <- sf::st_sfc(g)
    if (is.na(sf::st_crs(g)))
      stop("The sf object given as `x` has no CRS.", call. = FALSE)
    if (length(g) == 1L && all(sf::st_geometry_type(g) == "POINT")) {
      ll <- as.numeric(sf::st_coordinates(sf::st_transform(g, 4326)))[1:2]
      square(ll[1], ll[2])
    } else {
      box(as.numeric(sf::st_bbox(sf::st_transform(g, 4326))))
    }
  } else if (is.numeric(x) && length(x) %in% c(2L, 4L)) {
    lon <- x[c(TRUE, FALSE)]; lat <- x[c(FALSE, TRUE)]
    if (any(abs(lon) > 180) || any(abs(lat) > 90))
      stop(paste("Plain numbers in `x` are WGS84 longitude and latitude. For",
                 "projected coordinates, pass an sf point or bbox with its CRS."),
           call. = FALSE)
    if (length(x) == 2L) square(x[1], x[2])
    else box(c(min(lon), min(lat), max(lon), max(lat)))
  } else {
    stop(paste("`x` must be a longitude/latitude point c(lon, lat), a bbox",
               "c(xmin, ymin, xmax, ymax) in longitude/latitude, an sf point or",
               "geometry, or the path to a projected raster."), call. = FALSE)
  }

  if (a$extent[3] <= a$extent[1] || a$extent[4] <= a$extent[2])
    stop("The area given as `x` is empty.", call. = FALSE)
  ll <- gdalraster::bbox_transform(a$extent, a$crs, "EPSG:4326")
  lim <- 85.0511287798066
  if (ll[2] < -lim || ll[4] > lim)
    stop("Mapterhorn covers latitudes between 85.05 S and 85.05 N only.", call. = FALSE)
  list(extent = a$extent, crs = a$crs, ll = ll,
       res = if (isTRUE(a$grid)) a$res else NULL, grid = isTRUE(a$grid),
       point = isTRUE(a$point))
}

## XYZ tile index range covering a WGS84 bbox at zoom z.
.mt_tile_range <- function(ll, z) {
  n <- 2^z
  tx <- function(lon) min(max(floor((lon + 180) / 360 * n), 0), n - 1)
  ty <- function(lat) {
    r <- lat * pi / 180
    min(max(floor((1 - log(tan(r) + 1 / cos(r)) / pi) / 2 * n), 0), n - 1)
  }
  list(x = tx(ll[1]):tx(ll[3]), y = ty(ll[4]):ty(ll[2]))
}

.mt_tile_count <- function(ll, z) {
  r <- .mt_tile_range(ll, z)
  length(r$x) * length(r$y)
}

.mt_attribution <- function() {
  if (is.null(.mt_state$attribution)) {
    a <- tryCatch(jsonlite::fromJSON(.mt_attribution_url), error = function(e) e)
    if (inherits(a, "error"))
      stop(sprintf("Could not reach Mapterhorn at\n  %s\n%s", .mt_attribution_url,
                   conditionMessage(a)), call. = FALSE)
    .mt_state$attribution <- a[, c("source", "name", "producer", "license",
                                   "resolution", "website")]
  }
  .mt_state$attribution
}

## Sources whose coverage intersects the area, from Mapterhorn's coverage
## vector tiles, joined to their metadata. Read at the finest coverage zoom
## (<= 14) that needs no more than 16 tiles.
.mt_sources <- function(ll) {
  z <- 14L
  while (z > 0L && .mt_tile_count(ll, z) > 16L) z <- z - 1L
  r <- .mt_tile_range(ll, z)
  rect <- gdalraster::bbox_transform(ll, "EPSG:4326", "EPSG:3857")
  found <- character(0)
  fail <- function() stop(paste("Could not read Mapterhorn's coverage for this area.",
                                "The data is fetched on demand, so this needs an",
                                "internet connection."), call. = FALSE)
  for (tx in r$x) for (ty in r$y) {
    url <- sprintf(.mt_coverage_url, z, tx, ty)
    # open ocean has no coverage features at all, and the server answers 204;
    # tell that apart from a failed request, which must not silently drop a
    # source and hand back coarser data than exists
    status <- tryCatch(attr(curlGetHeaders(url), "status"), error = function(e) NA)
    if (identical(status, 204L)) next
    if (!identical(status, 200L)) fail()
    lyr <- tryCatch(methods::new(gdalraster::GDALVector, paste0("/vsicurl/", url),
                                 "coverage"), error = function(e) NULL)
    if (is.null(lyr)) fail()
    src <- tryCatch({
      lyr$setSpatialFilterRect(rect)
      unique(as.character(lyr$fetch(-1)$source))
    }, error = function(e) NULL, finally = lyr$close())
    if (is.null(src)) fail()
    found <- union(found, src)
  }
  att <- .mt_attribution()
  s <- att[att$source %in% found, , drop = FALSE]
  s <- s[order(s$resolution), , drop = FALSE]
  rownames(s) <- NULL
  s
}

## Zoom levels to mosaic, coarse to fine. Each source contributes the zoom it
## is built to at both latitude extremes of the area (they differ by one where
## the area crosses the threshold), capped at the zoom `res` needs; the global
## base underneath is always present.
.mt_levels <- function(sources, res, ll) {
  lat_near <- if (ll[2] <= 0 && ll[4] >= 0) 0 else min(abs(ll[c(2, 4)]))
  lat_far <- max(abs(ll[c(2, 4)]))
  need <- .mt_zoom_for(res, lat_near)
  z <- c(.mt_zoom_for(sources$resolution, lat_near),
         .mt_zoom_for(sources$resolution, lat_far))
  base <- min(.mt_planet_zoom, need)
  sort(unique(c(base, pmin(z[z > base], need))))
}

.mt_guard <- function(levels, ll, max_mb) {
  tiles <- sum(vapply(levels, function(z) .mt_tile_count(ll, z), numeric(1)))
  mb <- tiles * .mt_mb_per_tile
  if (mb <= max_mb) return(invisible(mb))
  stop(sprintf(paste("Reading this area would fetch about %s Mapterhorn tiles, an",
                     "estimated %.0f MB, above `max_mb = %s`. Reduce the area, use",
                     "a coarser `res`, or raise `max_mb` if you mean it."),
               format(tiles, big.mark = ","), mb, format(max_mb)), call. = FALSE)
}

## One zoom level as a Float32 elevation raster in EPSG:3857.
.mt_level_vrt <- function(z) {
  wms <- fs::file_temp(ext = "xml")
  writeLines(sprintf(paste0(
    '<GDAL_WMS><Service name="TMS"><ServerUrl>%s</ServerUrl></Service>',
    '<DataWindow><UpperLeftX>%.9f</UpperLeftX><UpperLeftY>%.9f</UpperLeftY>',
    '<LowerRightX>%.9f</LowerRightX><LowerRightY>%.9f</LowerRightY>',
    '<TileLevel>%d</TileLevel><TileCountX>1</TileCountX><TileCountY>1</TileCountY>',
    '<YOrigin>top</YOrigin></DataWindow><Projection>EPSG:3857</Projection>',
    '<BlockSizeX>512</BlockSizeX><BlockSizeY>512</BlockSizeY><BandsCount>3</BandsCount>',
    '<ZeroBlockHttpCodes>204,404</ZeroBlockHttpCodes></GDAL_WMS>'),
    .mt_tiles_url, -.mt_merc, .mt_merc, .mt_merc, -.mt_merc, z), wms)

  # Terrarium: elevation = R * 256 + G + B / 256 - 32768. A missing tile reads
  # as 0, 0, 0 and so decodes to exactly -32768, which becomes NoData.
  n <- 512 * 2^z
  px <- 2 * .mt_merc / n
  src <- function(b) sprintf(paste0('<SimpleSource><SourceFilename relativeToVRT="0">%s',
                                    '</SourceFilename><SourceBand>%d</SourceBand></SimpleSource>'),
                             wms, b)
  vrt <- fs::file_temp(ext = "vrt")
  writeLines(sprintf(paste0(
    '<VRTDataset rasterXSize="%.0f" rasterYSize="%.0f"><SRS>EPSG:3857</SRS>',
    '<GeoTransform>%.9f, %.12f, 0, %.9f, 0, %.12f</GeoTransform>',
    '<VRTRasterBand dataType="Float32" band="1" subClass="VRTDerivedRasterBand">',
    '<NoDataValue>-32768</NoDataValue><PixelFunctionType>expression</PixelFunctionType>',
    '<PixelFunctionArguments expression="B1 * 256 + B2 + B3 / 256 - 32768" dialect="muparser"/>',
    '%s%s%s</VRTRasterBand></VRTDataset>'),
    n, n, -.mt_merc, px, .mt_merc, -px, src(1), src(2), src(3)), vrt)
  vrt
}

## The lazy raster: each level warped onto the output grid on its own, then
## stacked coarse under fine.
##
## Resampling must never reach the encoded RGB. A VRT that resamples a derived
## band (buildVRT -r cubic, or any read at a buffer size other than the window)
## resamples the band's *sources* and decodes afterwards, and interpolating
## across a channel wrap - G every 1 m, R every 256 m - produces spikes that
## draw contour lines in a hillshade. Measured on 6 km at the Matterhorn:
## 40,053 spike cells from a cubic mosaic against 4,316 from the decoded tiles
## alone. The warper reads its source at native resolution, so warping each
## level directly means cubic only ever sees decoded heights; the warped levels
## share one grid, so stacking them resamples nothing.
.mt_vrt <- function(levels, extent, crs, res) {
  tr <- format(res, scientific = FALSE)
  lat <- mean(gdalraster::bbox_transform(extent, crs, "EPSG:4326")[c(2, 4)])
  warped <- vapply(levels, function(z) {
    out <- fs::file_temp(ext = "vrt")
    # Both of these keep reads window-invariant, so a metric's result cannot
    # depend on tile_size. Measured on a 1 km tile read whole and in 137 px
    # windows: without them values differed by up to 0.24 m.
    # * the warper derives the resampling kernel's scale from each request's
    #   own source and target window sizes, which differ slightly per window;
    #   fix it at the level's ground pixel (at mid-latitude) over res
    # * the default approximate transform is fitted per chunk; -et 0 is exact
    # -ovr NONE: the zoom is chosen here, never by the warper.
    scale <- sprintf("%.12f", .mt_px(z, lat) / res)
    gdalraster::warp(.mt_level_vrt(z), out, t_srs = crs,
      cl_arg = c("-of", "VRT", "-te", sprintf("%.6f", extent), "-tr", tr, tr,
                 "-r", "cubic", "-ot", "Float32", "-et", "0", "-ovr", "NONE",
                 "-wo", paste0("XSCALE=", scale), "-wo", paste0("YSCALE=", scale),
                 "-srcnodata", "-32768", "-dstnodata", "-9999"),
      quiet = TRUE)
    out
  }, character(1))
  if (length(warped) == 1L) return(warped)
  out <- fs::file_temp(ext = "vrt")
  gdalraster::buildVRT(out, warped,
    cl_arg = c("-srcnodata", "-9999", "-vrtnodata", "-9999"), quiet = TRUE)
  out
}

## Short, stable name for a CRS in cache file names. A hash of the WKT, not
## srs_find_epsg(), which can return a low-confidence wrong match.
.mt_crs_tag <- function(crs) {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(gdalraster::srs_to_wkt(crs), f)
  substr(unname(tools::md5sum(f)), 1, 8)
}

.mt_credit <- function(sources) {
  paste0("© Mapterhorn (", paste(sources$producer, collapse = "; "), ")")
}

## Resolve location, sources and resolution without touching any elevation data.
.mt_plan <- function(x, res, size, crs) {
  if (!is.null(res) && (!is.numeric(res) || length(res) != 1L || !is.finite(res) || res <= 0))
    stop("`res` must be a single positive number of metres, or NULL for the best available.",
         call. = FALSE)
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) || size <= 0)
    stop("`size` must be a single positive number of metres.", call. = FALSE)
  area <- .mt_area(x, size, crs)
  if (area$grid) {
    if (!is.null(res) && !isTRUE(all.equal(res, area$res)))
      stop("`res` can't differ from the resolution of the raster given as `x`.",
           call. = FALSE)
    res <- area$res
  }
  sources <- .mt_sources(area$ll)
  finest <- if (nrow(sources)) sources$resolution[1] else 30
  if (is.null(res)) {
    res <- finest
  } else if (res < finest * (1 - 1e-9)) {
    stop(sprintf(paste("The finest Mapterhorn data here is %s m (%s), so `res = %s`",
                       "would only interpolate. Use res >= %s."),
                 format(finest), sources$name[1], format(res), format(finest)),
         call. = FALSE)
  }
  extent <- if (area$grid) {
    area$extent
  } else if (area$point) {
    # keep a point's square exactly `size` wide: snap its corner, not both edges
    lo <- .data_snap(area$extent, res)[1:2]
    n <- ceiling(size / res - 1e-9) * res
    round(c(lo, lo + n), 6)
  } else {
    .data_snap(area$extent, res)
  }
  list(extent = extent, crs = area$crs, ll = area$ll, res = res, sources = sources)
}

#' Terrain for anywhere in the world, from Mapterhorn
#'
#' A digital terrain model for any place on Earth, from
#' [Mapterhorn](https://mapterhorn.com), which combines open elevation data into
#' one global dataset: 30 m Copernicus GLO-30 everywhere, and national LiDAR
#' terrain models where they have been published - 1 m across Germany, Austria
#' and Norway, 0.5 m in Switzerland, 0.25 m in parts of Baden-Württemberg, and
#' many more. Give a point to get a square around it, or an extent to get
#' exactly that area.
#'
#' Nothing is downloaded by default. The result is a small virtual raster that
#' reads Mapterhorn's tiles only when something reads it, and only the tiles
#' that read needs. Pass it to any function in this package.
#' `download = TRUE` stores a local copy instead, for offline work or repeated
#' sessions.
#'
#' @section Resolution:
#' `res = NULL` uses the finest source covering the area; see what that is with
#' [rvt_data_mapterhorn_sources()]. Asking for a finer `res` than that is an
#' error, because it could only interpolate. Where the area straddles sources
#' of different resolution - the edge of a national LiDAR survey, say - the
#' coarser source fills in wherever the finer one has no data, so the result
#' has no holes.
#'
#' Heights are stored to 1/256 m (about 4 mm), which is below the accuracy of
#' any source.
#'
#' @section Coordinate system:
#' Mapterhorn's tiles are in Web Mercator, whose map units are stretched by
#' 1/cos(latitude): 1.6 times at 51 degrees north. Since every distance in this
#' package is in map units, the data is always reprojected to a coordinate
#' system in metres - by default the UTM zone of the area's centre (standard
#' zones, e.g. EPSG:32632), or whatever `crs` you pass. Pass the national
#' system to line data up with other sources, e.g. `crs = 25832` for Germany.
#'
#' @section Location:
#' * **A point** returns a `size` x `size` m square centred on it. Plain
#'   numbers `c(lon, lat)` are WGS84; an sf point uses its own coordinate system.
#' * **An extent** returns exactly that area, snapped outward to whole cells:
#'   an `sf::st_bbox()`, any sf geometry, or plain WGS84
#'   `c(xmin, ymin, xmax, ymax)`.
#' * **A raster path** returns data on exactly that raster's grid and
#'   coordinate system, which must be projected.
#'
#' @section Licence:
#' Mapterhorn's data comes from many producers under their own open licences
#' (CC BY 4.0, Datenlizenz Deutschland, and others). The sources used are
#' listed in a message the first time each is fetched in a session, and by
#' [rvt_data_mapterhorn_sources()]. Credit Mapterhorn and each source; see
#' <https://mapterhorn.com/attribution> for the full terms.
#'
#' @param x location: a WGS84 `c(lon, lat)` point or `c(xmin, ymin, xmax,
#'   ymax)` extent, an sf point, geometry or bbox, or the path to a raster in a
#'   projected coordinate system
#' @param res output resolution in metres, or `NULL` (default) for the finest
#'   source available
#' @param size side of the square returned for a point, in metres (default 1000)
#' @param crs output coordinate system: an EPSG code, `"EPSG:..."` string or
#'   WKT, in metres. `NULL` (default) uses the local UTM zone.
#' @param download `FALSE` (default) returns a virtual raster reading the
#'   tiles on demand; `TRUE` stores a local copy in the user cache
#' @param max_mb stop before reading if the estimated download exceeds this
#'   many megabytes (default 500)
#' @param refresh with `download = TRUE`, fetch again even if cached
#' @inheritParams rvt_svf
#' @return path to the raster (a VRT, or a Cloud-Optimized GeoTIFF with
#'   `download = TRUE`), so it pipes into any function
#' @seealso [rvt_data_mapterhorn_sources()], and [rvt_data_lgln()] for Lower
#'   Saxony's terrain, surface model and orthophotos at their native grid
#' @examples
#' \donttest{
#' # the Matterhorn, 3 km across at Switzerland's 0.5 m, shown at 2 m
#' rvt_data_mapterhorn(c(7.6586, 45.9763), res = 2, size = 3000) |>
#'   rvt_hillshade() |>
#'   rvt_plot()
#'
#' rvt_data_mapterhorn_sources(c(7.6586, 45.9763))
#' }
#' @export
rvt_data_mapterhorn <- function(x, res = NULL, size = 1000, crs = NULL,
                                download = FALSE, max_mb = 500, refresh = FALSE,
                                threads = rvt_threads()) {
  p <- .mt_plan(x, res, size, crs)
  levels <- .mt_levels(p$sources, p$res, p$ll)
  .mt_guard(levels, p$ll, max_mb)

  if (isTRUE(download)) {
    out_path <- .as_path(fs::path(tools::R_user_dir("rvtr", "cache"), "mapterhorn",
      sprintf("mapterhorn_%s_%sm_%s.tif", .mt_crs_tag(p$crs),
              format(p$res, scientific = FALSE),
              paste(format(round(p$extent, 6), scientific = FALSE), collapse = "_"))))
    if (!refresh && fs::file_exists(out_path)) return(out_path)
  }

  new <- p$sources[!p$sources$source %in% .mt_state$credited, , drop = FALSE]
  if (nrow(new)) {
    message("Mapterhorn terrain, https://mapterhorn.com/attribution. Sources here:\n",
            paste0("  ", new$name, " - ", new$producer, ", ", new$license,
                   collapse = "\n"))
    .mt_state$credited <- union(.mt_state$credited, new$source)
  }

  vrt <- tryCatch(.mt_vrt(levels, p$extent, p$crs, p$res),
                  error = function(e) stop(paste("Could not open the Mapterhorn tiles:\n",
                                                 conditionMessage(e)), call. = FALSE))
  if (!isTRUE(download)) return(.as_path(vrt))

  md <- c(SOURCE = paste("Mapterhorn:", paste(p$sources$name, collapse = "; ")),
          ATTRIBUTION = .mt_credit(p$sources),
          LICENSE = paste(unique(p$sources$license), collapse = "; "))
  if (!.data_write_cog(vrt, out_path, "FLOATING_POINT", md, threads))
    stop(paste("Could not download the Mapterhorn terrain. The data is fetched on",
               "demand, so this needs an internet connection."), call. = FALSE)
  out_path
}

#' Which Mapterhorn sources cover a place
#'
#' Lists the elevation sources Mapterhorn holds for an area, finest first,
#' without fetching any elevation data - to see what resolution
#' [rvt_data_mapterhorn()] will return, and whom to credit.
#'
#' @inheritParams rvt_data_mapterhorn
#' @return a data frame with columns `source`, `name`, `producer`, `license`,
#'   `resolution` (metres) and `website`
#' @seealso [rvt_data_mapterhorn()]
#' @examples
#' \donttest{
#' rvt_data_mapterhorn_sources(c(9.9464, 51.6317))
#' }
#' @export
rvt_data_mapterhorn_sources <- function(x, size = 1000) {
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) || size <= 0)
    stop("`size` must be a single positive number of metres.", call. = FALSE)
  .mt_sources(.mt_area(x, size, NULL)$ll)
}
