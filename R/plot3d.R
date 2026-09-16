## Interactive 3D view: the terrain as a surface, a finished raster draped
## over it as a texture.
##
## The relief comes from the *drape*, not from rgl's lighting - this package
## already computes better shading than a single OpenGL lamp, so the surface is
## deliberately unlit (lit = FALSE) and shows exactly the image handed to it.
##
## Everything stays in the raster's own coordinate system. Map units are metres
## for the projected CRSs this package works with, so x, y and z share units
## and vertical exaggeration is one multiplication - no reprojection anywhere.

## The DEM read downsampled to max_dim, plus the cell centres in map units.
.plot3d_grid <- function(path, max_dim, band = 1L) {
  ds <- methods::new(gdalraster::GDALRaster, path, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  scale <- min(1, max_dim / max(nx, ny))
  out_nx <- max(2L, round(nx * scale))
  out_ny <- max(2L, round(ny * scale))

  # NoData comes back as NA already, whether the band flags NaN (the bundled
  # tile) or -9999 (everything this package writes), so holes become gaps in
  # the surface with no masking of our own
  z <- gdalraster::read_ds(ds, bands = band, out_xsize = out_nx,
                           out_ysize = out_ny)

  bb <- ds$bbox()
  mid <- function(v) (v[-length(v)] + v[-1]) / 2
  x <- mid(seq(bb[1], bb[3], length.out = out_nx + 1L))
  y <- mid(seq(bb[2], bb[4], length.out = out_ny + 1L))

  # read_ds() returns values row by row, north row first, west to east within
  # a row. matrix() fills down columns, so filling with nrow = out_nx puts a
  # raster *row* into each column: the result is indexed [y, x], the transpose
  # of what surface3d() wants. Hence fill with one column per raster row
  # (nrow = out_nx is wrong; use out_nx as the row count only after t()).
  #
  # Verified against the raster, not assumed: read_ds()[1:5] is the north row
  # west to east. Getting this wrong does not mirror the map, it tears the
  # surface into folds, and a square test grid hides it completely.
  zm <- t(matrix(z, nrow = out_ny, ncol = out_nx, byrow = TRUE))  # [x, y], y N->S
  zm <- zm[, rev(seq_len(out_ny)), drop = FALSE]                  # y ascending

  list(x = x, y = y, z = zm, nx = out_nx, ny = out_ny, bbox = bb)
}

## The drape as a PNG the size of the source raster's own pixels, so the
## texture is never coarser than what the user asked to see.
.plot3d_texture <- function(drape, col, max_dim) {
  if (inherits(drape, "rvt_stack")) drape <- rvt_render(drape)
  drape <- .as_path(drape)
  ds <- methods::new(gdalraster::GDALRaster, drape, read_only = TRUE)
  on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize()
  scale <- min(1, max_dim / max(nx, ny))
  out_nx <- max(2L, round(nx * scale))
  out_ny <- max(2L, round(ny * scale))
  nb <- ds$getRasterCount()
  bands <- if (nb >= 3L) 1:3 else 1L

  v <- gdalraster::read_ds(ds, bands = bands, out_xsize = out_nx,
                           out_ysize = out_ny)
  png <- fs::file_temp(ext = "png")

  if (length(bands) == 3L) {
    # an orthophoto or a rendered blend: already colours, only the 0-255 scale
    # to undo
    a <- array(pmin(pmax(as.numeric(v), 0), 255) / 255, c(out_nx, out_ny, 3L))
    rgb <- grDevices::rgb(a[, , 1], a[, , 2], a[, , 3])
    rgb[is.na(as.numeric(v[seq_len(out_nx * out_ny)]))] <- "#FFFFFF"
    m <- matrix(rgb, out_nx, out_ny)
  } else {
    # a single-band metric: stretch to its own range and map through the
    # palette, the same convention rvt_plot() uses
    rng <- range(v, na.rm = TRUE)
    # .resolve_col() hands back the colours themselves, not a ramp function
    pal <- .resolve_col(col, 256L)
    idx <- if (diff(rng) > 0) {
      1L + round((v - rng[1]) / diff(rng) * 255)
    } else {
      rep(128L, length(v))
    }
    m <- matrix(pal[idx], out_nx, out_ny)
    m[is.na(m)] <- "#FFFFFF"
  }

  # A plain grDevices::png() writes a *paletted* PNG when the image has few
  # enough colours, which makes a poor texture; ragg always writes RGB. The
  # matrix is [x, y] with y running north to south, and as.raster() takes rows
  # top to bottom, so it transposes to put north at the top.
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(png, width = out_nx, height = out_ny, background = "white")
  } else {
    grDevices::png(png, width = out_nx, height = out_ny, type = "cairo")
  }
  op <- graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  graphics::plot.new()
  graphics::plot.window(c(0, 1), c(0, 1))
  graphics::rasterImage(grDevices::as.raster(t(m)), 0, 0, 1, 1,
                        interpolate = FALSE)
  graphics::par(op)
  grDevices::dev.off()
  png
}

#' Interactive 3D view of a terrain model
#'
#' Draws the terrain as a surface you can rotate and zoom, with a finished
#' raster draped over it: a blend from [rvt_blend()], a relit orthophoto, an
#' orthophoto, or any single metric coloured by a palette. Needs the `rgl`
#' package.
#'
#' The shading you see comes from the **drape**, not from a light in the 3D
#' scene. This package computes relief better than a single OpenGL lamp can, so
#' the surface is drawn unlit and shows exactly the image given to it. Build
#' the picture first, in 2D, then drape it.
#'
#' Everything stays in the raster's own coordinate system, which is not
#' reprojected. For the projected systems this package is built for, map units
#' are metres, so heights and distances are already comparable and
#' `exaggeration` is a plain multiplier on the heights.
#'
#' ```r
#' # the terrain itself, coloured by height: col is the second argument,
#' # exactly as in rvt_plot()
#' dem |> rvt_plot3d("Terrain 2")
#'
#' # an orthophoto, a relit image, or any other finished raster
#' rvt_plot3d(dem, drape = ortho)
#'
#' # for a page or a vignette rather than a window
#' rvt_plot3d(dem, drape = ortho, widget = TRUE)
#' ```
#'
#' @section Size:
#' The DEM is read downsampled to `max_dim`, exactly as [rvt_plot()] does, so a
#' huge raster previews cheaply. One cell becomes one vertex, so the default
#' 1000 is already a million of them - enough for a 1 km tile at 1 m. Raising
#' it much further will make navigation stutter long before it runs out of
#' memory, and `widget = TRUE` pages grow quickly, since the vertices are
#' written into the HTML.
#'
#' @param dem path to a terrain or surface model, or a vector of paths to
#'   mosaic
#' @param col palette for the terrain, or for a single-band `drape`. Second
#'   argument, as in [rvt_plot()], so `dem |> rvt_plot3d("Terrain 2")` works
#'   the same way in both. [rvt_palettes()] shows them all
#' @param drape what to lay over the surface: a raster path, an `rvt_stack`
#'   from [rvt_blend()], or `NULL` (default) to colour the terrain itself by
#'   `col`. Three-band rasters are drawn as colours; single-band ones are
#'   stretched and put through the palette
#' @param max_dim longest side, in pixels, to read at (default 1000)
#' @param exaggeration multiplier on the heights (default 1, true shape)
#' @param background background colour of the scene
#' @param widget return an `rglwidget()` for an HTML page instead of opening a
#'   window (default `FALSE`)
#' @return `dem`, invisibly, so it pipes; or the widget when `widget = TRUE`
#' @seealso [rvt_plot()] for the 2D version, [rvt_blend()] to build a drape
#' @examples
#' \dontrun{
#' dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
#'
#' # a composite draped over the terrain, with the relief exaggerated
#' stack <- rvt_hillshade(dem) |>
#'   rvt_blend(rvt_svf(dem), "multiply", opacity = 0.25)
#' rvt_plot3d(dem, drape = stack, exaggeration = 2)
#' }
#' @export
rvt_plot3d <- function(dem, col = "Grays", drape = NULL, max_dim = 1000,
                        exaggeration = 1, background = "white",
                        widget = FALSE) {
  if (!requireNamespace("rgl", quietly = TRUE))
    stop("rvt_plot3d() needs the rgl package. Install it with install.packages(\"rgl\").",
         call. = FALSE)
  if (!is.numeric(exaggeration) || length(exaggeration) != 1L ||
      !is.finite(exaggeration) || exaggeration <= 0)
    stop("`exaggeration` must be a single positive number.", call. = FALSE)
  if (!is.numeric(max_dim) || length(max_dim) != 1L || max_dim < 2)
    stop("`max_dim` must be a single number of at least 2.", call. = FALSE)

  # `col` is second, matching rvt_plot(), so a raster passed positionally is
  # the easy slip in this direction - say so rather than letting it fail
  # inside the palette lookup
  if (inherits(col, "rvt_stack") ||
      (is.character(col) && length(col) == 1L && fs::file_exists(.as_path(col)))) {
    stop(paste("`col` is a palette, as in rvt_plot(). To lay a raster over the",
               "terrain, name it: rvt_plot3d(dem, drape = ...)."), call. = FALSE)
  }

  dem_path <- rvt_mosaic(dem)
  g <- .plot3d_grid(dem_path, max_dim)

  # with no drape the terrain colours itself, which needs no texture file
  texture <- if (is.null(drape)) NULL else .plot3d_texture(drape, col, max_dim)
  on.exit(if (!is.null(texture)) .rm_path(texture), add = TRUE)

  rgl::open3d()
  rgl::bg3d(color = background)
  if (is.null(texture)) {
    pal <- .resolve_col(col, 256L)
    rng <- range(g$z, na.rm = TRUE)
    idx <- if (diff(rng) > 0) 1L + round((g$z - rng[1]) / diff(rng) * 255) else
      array(128L, dim(g$z))
    rgl::surface3d(g$x, g$y, g$z * exaggeration, color = pal[idx],
                   lit = FALSE, back = "filled")
  } else {
    # Texture coordinates run 0-1 across the grid: s west to east, and t
    # **0 at the south edge, 1 at the north**. rgl maps t = 1 to the PNG's
    # first row, which is the north edge, so t must ascend with y - the
    # opposite way round shipped once and drew the whole drape upside down.
    # Measured, not reasoned: against the same surface coloured per vertex,
    # t ascending correlates 0.988 and t descending -0.784.
    s_mat <- matrix(seq(0, 1, length.out = g$nx), g$nx, g$ny)
    t_mat <- matrix(rep(seq(0, 1, length.out = g$ny), each = g$nx), g$nx, g$ny)
    rgl::surface3d(g$x, g$y, g$z * exaggeration, texture = texture,
                   texture_s = s_mat, texture_t = t_mat, color = "white",
                   lit = FALSE, back = "filled")
  }

  # Heights are true metres, the same units as x and y, so the *data* decides
  # the shape: aspect3d(1, 1, 1) would force a cube and squash a 1 km tile
  # with 70 m of relief into a near-flat plate - which is exactly what made
  # early renders unreadable. Passing the spans keeps the terrain in
  # proportion, and `exaggeration` stays the only vertical control.
  rgl::aspect3d(diff(range(g$x)), diff(range(g$y)),
                diff(range(g$z * exaggeration, na.rm = TRUE)))

  if (isTRUE(widget)) return(rgl::rglwidget())
  invisible(dem)
}
