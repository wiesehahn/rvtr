# Interactive 3D view of a terrain model

Draws the terrain as a surface you can rotate and zoom, with a finished
raster draped over it: a blend from
[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md),
a relit orthophoto, an orthophoto, or any single metric coloured by a
palette. Needs the `rgl` package.

## Usage

``` r
rvt_plot3d(
  dem,
  col = "Grays",
  drape = NULL,
  max_dim = 1000,
  exaggeration = 1,
  background = "white",
  widget = FALSE
)
```

## Arguments

- dem:

  path to a terrain or surface model, or a vector of paths to mosaic

- col:

  palette for the terrain, or for a single-band `drape`. Second
  argument, as in
  [`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md),
  so `dem |> rvt_plot3d("Terrain 2")` works the same way in both.
  [`rvt_palettes()`](https://wiesehahn.github.io/rvtr/reference/rvt_palettes.md)
  shows them all

- drape:

  what to lay over the surface: a raster path, an `rvt_stack` from
  [`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md),
  or `NULL` (default) to colour the terrain itself by `col`. Three-band
  rasters are drawn as colours; single-band ones are stretched and put
  through the palette

- max_dim:

  longest side, in pixels, to read at (default 1000)

- exaggeration:

  multiplier on the heights (default 1, true shape)

- background:

  background colour of the scene

- widget:

  return an `rglwidget()` for an HTML page instead of opening a window
  (default `FALSE`)

## Value

`dem`, invisibly, so it pipes; or the widget when `widget = TRUE`

## Details

The shading you see comes from the **drape**, not from a light in the 3D
scene. This package computes relief better than a single OpenGL lamp
can, so the surface is drawn unlit and shows exactly the image given to
it. Build the picture first, in 2D, then drape it.

Everything stays in the raster's own coordinate system, which is not
reprojected. For the projected systems this package is built for, map
units are metres, so heights and distances are already comparable and
`exaggeration` is a plain multiplier on the heights.

    # the terrain itself, coloured by height: col is the second argument,
    # exactly as in rvt_plot()
    dem |> rvt_plot3d("Terrain 2")

    # an orthophoto, a relit image, or any other finished raster
    dem |> rvt_plot3d(drape = ortho)

    # for a page or a vignette rather than a window
    dem |> rvt_plot3d(drape = ortho, widget = TRUE)

## Size

The DEM is read downsampled to `max_dim`, exactly as
[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)
does, so a huge raster previews cheaply. One cell becomes one vertex, so
the default 1000 is already a million of them - enough for a 1 km tile
at 1 m. Raising it much further will make navigation stutter long before
it runs out of memory, and `widget = TRUE` pages grow quickly, since the
vertices are written into the HTML.

## See also

[`rvt_plot()`](https://wiesehahn.github.io/rvtr/reference/rvt_plot.md)
for the 2D version,
[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md)
to build a drape

## Examples

``` r
if (FALSE) { # \dontrun{
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")

# a composite draped over the terrain, with the relief exaggerated
stack <- dem |> rvt_hillshade() |>
  rvt_blend(rvt_svf(dem), "multiply", opacity = 0.25)
dem |> rvt_plot3d(drape = stack, exaggeration = 2)
} # }
```
