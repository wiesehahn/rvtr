# Figure conventions shared by every vignette. Sourced from each vignette's
# hidden options chunk; knitr runs in vignettes/, so the relative path works
# for pkgdown and R CMD build alike.
#
# Figures are sized in pixels. Every image panel is 400 x 400 px, with a 36 px
# title strip above and a 2 px gutter around it, so a panel cell is 404 x 440 px
# and a 1 x 3 figure is 1212 x 440 px. Rasters drawn with max_dim = 400 are
# therefore shown at their own pixels, and panels match across vignettes.
#
# dpi does not change the pixels, only how large text is drawn within them. It
# was chosen by rendering the same 1 x 3 figure at 100, 120 and 150 dpi and
# viewing each at a typical 730 px page width: the files were the same size
# (263-267 kB), but titles were unreadable at 100 and small at 120. 150 it is.

vignette_dpi <- 150
panel_px <- 400
title_px <- 36
gutter_px <- 2

# pixels -> the inches knitr's fig.width / fig.height expect
# (a quarter pixel added: 440 / 150 in is 2.9333..., and the device rounds
# width x res down, which would lose the last row)
px <- function(n) (n + 0.25) / vignette_dpi

# figure size for a grid of panels, in inches for the chunk header. Grids with
# more than three columns pass a smaller `panel` (e.g. 240 for five), so no
# figure is wider than a 1 x 3 row: the page shows ~730 px, and a 2020 px wide
# 3 x 5 grid of 400 px panels cost 1.46 MB for detail nobody could see. Draw
# those panels with the same max_dim.
fig_w <- function(ncol, panel = panel_px) px(ncol * (panel + 2 * gutter_px))
fig_h <- function(nrow, panel = panel_px) px(nrow * (panel + title_px + 2 * gutter_px))

# JPEG through ragg: sharp text, quality 80 keeps rasters and photos small
# without visible artefacts. knitr knows ragg_png by name but not a ragg JPEG
# device, hence this wrapper.
ragg_jpeg <- function(filename, width, height, ...) {
  ragg::agg_jpeg(filename, width = width, height = height, units = "in",
                 res = vignette_dpi, quality = 80)
}

knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", out.width = "100%",
  fig.width = fig_w(3), fig.height = fig_h(1), dpi = vignette_dpi,
  dev = "ragg_jpeg", fig.ext = "jpeg",
  panels = TRUE
)

# Panel layout comes from the chunk option `layout = c(rows, cols)`, not a
# par(mfrow = ...) in the visible code. Setting mfrow makes R shrink all text
# to 83% for a 2 x 2 grid and to 66% as soon as there are three rows or
# columns, so titles came out smaller in every 1 x 3 figure than in a 1 x 2
# one. Here `cex = 1` follows mfrow in the same call and undoes that, so every
# title is the same size. Margins are in pixels, not text lines, so each
# panel's image area is exactly panel_px square whatever the dpi.
knitr::knit_hooks$set(panels = function(before, options, envir) {
  if (before) graphics::par(
    mfrow = if (is.null(options$layout)) c(1, 1) else options$layout,
    cex = 1,
    mai = c(gutter_px, gutter_px, title_px + gutter_px, gutter_px) / vignette_dpi,
    font.main = 1, cex.main = 1, col.main = "grey20")
})
