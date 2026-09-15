## A browser for the colour palettes rvt_plot() accepts by name.
##
## Base R ships 115 HCL palettes but no function to look at them - only a
## swatch snippet inside the examples of ?hcl.colors. colorspace has one, but
## is not a dependency, and a palette is easier chosen by eye than by name.

.palette_types <- c("sequential", "diverging", "divergingx", "qualitative")

## What rvt_plot()'s `col` means: a single string is a palette name, anything
## longer is colours as given. Resolved here rather than by handing the name to
## hcl.colors(), whose errors are translated into the session language and
## which rejects a plain "red" as an ambiguous palette prefix without saying
## which palettes it could have meant.
.resolve_col <- function(col, n = 256L) {
  if (is.function(col)) return(col)
  if (!is.character(col) || length(col) == 0L)
    stop("`col` must be colours or the name of a palette - see rvt_palettes().",
         call. = FALSE)

  if (length(col) > 1L) {
    ok <- tryCatch({ grDevices::col2rgb(col); TRUE }, error = function(e) FALSE)
    if (!ok) stop("`col` contains something that is not a colour.", call. = FALSE)
    return(col)
  }

  pals <- grDevices::hcl.pals()
  # case, spaces and hyphens ignored, as hcl.colors() itself does
  key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  hit <- pals[key(pals) == key(col)]
  if (length(hit) == 0L) hit <- pals[startsWith(key(pals), key(col))]

  if (length(hit) == 1L) return(grDevices::hcl.colors(n, hit))
  if (length(hit) > 1L)
    stop(sprintf('`col = "%s"` matches several palettes: %s. See rvt_palettes().',
                 col, paste(hit, collapse = ", ")), call. = FALSE)
  stop(sprintf('`col = "%s"` is not a palette name. See rvt_palettes() for the %d available.',
               col, length(pals)), call. = FALSE)
}

#' Browse the colour palettes
#'
#' Draws every palette [rvt_plot()] accepts by name as a colour ramp, grouped
#' by kind, so you can pick one by eye rather than by guessing a name. Pass the
#' name you like straight to `rvt_plot()`: `rvt_plot(svf, "Viridis")`.
#'
#' Which kind to reach for depends on what the values mean:
#'
#' * **sequential** for magnitudes running one way - sky-view factor, slope,
#'   openness, elevation. Most relief visualizations want one of these, and a
#'   plain `"Grays"` is hard to beat for reading shape.
#' * **diverging** for signed results centred on zero - curvature, [rvt_log()],
#'   [rvt_slrm()], [rvt_msrm()], [rvt_dev()] - so convex and concave read as
#'   opposites. Keep the stretch symmetric about zero, as [rvt_plot()]
#'   explains, or zero will not land on the neutral middle.
#' * **divergingx** join two sequential ramps with no requirement that they
#'   balance. Many are familiar from other tools (`"RdBu"`, `"PuOr"`,
#'   `"Spectral"`), but one arm can end much darker than the other, so equal
#'   positive and negative values need not look equally strong: `"PuOr"`
#'   ends at lightness 10 on its purple side against 31 on its orange side.
#'   `"Cividis"` is filed here too, although it runs dark to light and is
#'   really sequential. Prefer **diverging** when that balance matters.
#' * **qualitative** for classes with no order, such as [rvt_geomorphons()]:
#'   a ramp would imply an order the classes do not have.
#'
#' @param type `NULL` (default) for all of them, or one of `"sequential"`,
#'   `"diverging"`, `"divergingx"`, `"qualitative"`
#' @param n number of colours drawn in each ramp (default 64)
#' @return the palette names shown, invisibly
#' @seealso [rvt_plot()]
#' @examples
#' rvt_palettes("diverging")
#' names <- rvt_palettes("sequential")
#' @export
rvt_palettes <- function(type = NULL, n = 64) {
  types <- .palette_types
  if (!is.null(type)) {
    if (!is.character(type) || length(type) != 1L || !type %in% .palette_types)
      stop(paste('`type` must be NULL or one of "sequential", "diverging",',
                 '"divergingx", "qualitative".'), call. = FALSE)
    types <- type
  }
  if (!is.numeric(n) || length(n) != 1L || n < 2)
    stop("`n` must be a single number of at least 2", call. = FALSE)
  n <- as.integer(n)

  # one row per palette, plus a heading row per kind when showing several
  rows <- list()
  for (t in types) {
    if (length(types) > 1L) rows[[length(rows) + 1L]] <- list(header = t)
    for (nm in grDevices::hcl.pals(t)) rows[[length(rows) + 1L]] <- list(name = nm)
  }

  # at most 30 rows a column, but no more rows than there are, so a single
  # kind fills the height instead of leaving the bottom of the device empty
  per_col <- min(30L, length(rows))
  ncol <- ceiling(length(rows) / per_col)
  op <- graphics::par(mfrow = c(1L, ncol), mar = rep(0.2, 4), cex = 0.7)
  on.exit(graphics::par(op), add = TRUE)

  # names get the left third, right-aligned against the ramp; the ramp the rest
  name_edge <- 0.34
  x <- seq(0.36, 1, length.out = n + 1L)
  for (k in seq_len(ncol)) {
    idx <- ((k - 1L) * per_col + 1L):min(k * per_col, length(rows))
    graphics::plot.new()
    graphics::plot.window(c(0, 1), c(0, per_col))
    for (j in seq_along(idx)) {
      r <- rows[[idx[j]]]
      y <- per_col - j                 # top-down
      if (!is.null(r$header)) {
        graphics::text(x[1L], y + 0.35, toupper(r$header), adj = c(0, 0.5), font = 2)
      } else {
        graphics::rect(x[-(n + 1L)], y, x[-1L], y + 0.7,
                       col = grDevices::hcl.colors(n, r$name), border = NA)
        graphics::text(name_edge, y + 0.35, r$name, adj = c(1, 0.5))
      }
    }
  }
  invisible(unlist(lapply(types, grDevices::hcl.pals), use.names = FALSE))
}
