# rvtr

Topographic visualizations for R, following the Relief Visualization Toolbox
(Kokalj & Somrak 2019), but with no Python dependency and no WhiteboxTools:
sky-view factor, positive and negative openness, local dominance, the
multi-scale relief model (MSRM), and the archaeological **VAT** blend.

Built for large inputs: a fused, OpenMP-parallel C++ horizon kernel, tile-by-tile
processing so memory stays bounded, and mosaics of many rasters handled through a
GDAL virtual raster so tiles never produce seams.

> The package is named `rvtr` rather than `relief_vizualization_toolbox` because R
> package names cannot contain underscores. The repository directory keeps the
> longer name.

## Install

Open `relief_vizualization_toolbox.Rproj` in RStudio, then in the console:

```r
install.packages(".", repos = NULL, type = "source")
```

Needs a C++ compiler to build (on Windows, install Rtools first — RStudio
prompts for this automatically if it's missing). Everything else, including
GDAL via `gdalraster`, installs as a regular dependency.

## Use

Runs as-is against the 1 m sample tile bundled with the package. Output path is
optional — omit it and you get a temp file; every function returns the path it
wrote:

```r
library(rvtr)

dem <- system.file("extdata", "dtm1.tif", package = "rvtr")

rvt_hillshade(dem)           # shaded relief, 0-1
rvt_multi_hillshade(dem)     # shaded relief lit from four sides at once
rvt_shadow(dem)              # cast shadow, 1 lit / 0 shadowed
rvt_daylight(dem)            # hours of direct sun per day, over a year
rvt_insolation(dem)          # clear-sky solar energy, kWh/m2 per year
rvt_slope(dem)               # steepness, degrees
rvt_svf(dem)                 # sky-view factor, 0-1
rvt_asvf(dem)                # anisotropic sky-view factor, 0-1
rvt_sky_illumination(dem)    # diffuse sky illumination, flat ground = 1
rvt_openness(dem)            # positive openness, degrees
rvt_openness_negative(dem)   # negative openness, degrees
rvt_local_dominance(dem)     # local dominance
rvt_slrm(dem)                # simple local relief model, metres (= rvt_tpi)
rvt_dev(dem)                 # deviation from mean elevation, standardised
rvt_msrm(dem)                # multi-scale relief model, metres
rvt_mstp(dem)                # multi-scale topographic position, RGB
rvt_vat(dem)                 # archaeological VAT (combined), 0-1

rvt_aspect(dem)              # direction each cell faces, degrees
rvt_curvature(dem)           # profile curvature; convex +, concave -
rvt_log(dem, sigma = 3)      # Laplacian of Gaussian, edges at a chosen scale
rvt_geomorphons(dem)         # landform class per cell (10 classes)
```

For your own data, swap `dem` for a path to your DEM, and pass an output path
if you want to keep the file, e.g. `rvt_svf(dem, "svf.tif")`. The DEM is always
the first argument, so these also work with the pipe: `dem |> rvt_vat("vat.tif")`.

**Every distance is in map units** — metres, normally — never pixels. So
`rvt_svf(dem, reach = 10)` looks 10 m out whether the DEM is 2 m or 0.25 m, and
the radii published in the literature can be typed in as they stand. Distances
snap to whole cells: a request smaller than one cell is an error rather than a
silent round-up, and anything the grid can't hit closely says so.

### Which one to reach for

| function | what it shows | flat ground reads |
|---|---|---|
| `rvt_hillshade()` | classic shaded relief from one sun position; intuitive, but features along the light fade out | varies |
| `rvt_multi_hillshade()` | the same lit from four directions at once — nothing is hidden and no side is left in deep shadow, at the cost of some contrast | varies |
| `rvt_shadow()` | true cast shadow for a sun position — unlike hillshade, terrain in the way actually blocks the light | 1 (lit) |
| `rvt_daylight()` | hours of direct sun per day, integrated over the year (or any day range) — cast shadow and slope aspect together | day length |
| `rvt_insolation()` | clear-sky solar **energy** in kWh/m², beam and diffuse, following GRASS `r.sun`'s ESRA model — what the hour counts are a proxy for | varies |
| `rvt_slope()` | steepness alone; very sharp on breaks of slope, and noisy on noisy DEMs | 0° |
| `rvt_svf()` | how much sky each cell sees — like an overcast sky with no sun direction, so nothing hides in shadow | 1 |
| `rvt_asvf()` | sky-view factor under a sky brighter on one side: directional like a hillshade, but nothing is lost to shadow | 1 |
| `rvt_sky_illumination()` | diffuse daylight actually received, counting both the horizon and which way the ground tilts | 1 |
| `rvt_openness()` | the same horizon on a wider scale, spreading convex and concave apart; banks and mounds bright | 90° |
| `rvt_openness_negative()` | the same on inverted terrain, so ditches, pits and hollow ways come forward instead | 90° |
| `rvt_local_dominance()` | how far the surroundings fall away below an observer standing there; good at low, spread-out features | 1 |
| `rvt_slrm()` | height above/below the smoothed terrain **in metres**, at one scale you choose | 0 |
| `rvt_msrm()` | the same idea across a range of scales at once, so you needn't know the feature size up front | 0 |
| `rvt_dev()` | the same standardised by local roughness, so a bank on flat ground and a terrace on a rough hillside read alike | 0 |
| `rvt_mstp()` | RGB composite colouring each cell by the *scale* at which it stands out: blue fine, green mid, red broad | dark |
| `rvt_curvature()` | how sharply the ground bends — the natural detector for breaks of slope, but noisy | 0 |
| `rvt_log()` | the same idea with a scale knob: edges at the size you ask for, noise below it suppressed | 0 |
| `rvt_geomorphons()` | a landform **class** per cell (ridge, spur, hollow, valley…) rather than a number | flat |
| `rvt_vat()` | ready-made greyscale composite of four of the above, tuned for spotting earthworks | — |

Each help page (`?rvt_svf` and so on) has a **How it works** section covering the
sampling geometry, a **Tuning** section on adapting the parameters, and a
**Recommended settings** section drawing on Kokalj & Hesse's *Guide to Good
Practice* — what each metric suits, what it can't show, and the values and
display stretches that work in flat, moderate and steep terrain. Their radii
are in metres, and so is every distance these functions take, so the published
values can be used as they stand.

`?rvtr` collects the guide's advice on **which** visualization to reach for in
which terrain, and in what order. The short version: start with a hillshade
for orientation, expect to need several techniques rather than one, and pick
by terrain — trend removal and local dominance on flat ground, sky-view factor
once there are real slopes, openness where the topography is steep enough that
other methods saturate.

> Kokalj, Ž. & Hesse, R. (2017) *Airborne Laser Scanning Raster Data
> Visualization: A Guide to Good Practice*. Prostor, kraj, čas 14. Ljubljana:
> Založba ZRC. <https://doi.org/10.3986/9789612549848> (open access)

### Mosaics

Pass a vector of paths. They are wrapped in a VRT, so the horizon search reads
across file boundaries and the result is identical to processing the same area as
one raster (verified bit-for-bit against a single-file run):

```r
tiles <- list.files("dtm_tiles", "\\.tif$", full.names = TRUE)
tiles |> rvt_vat("vat_mosaic.tif")
```

### Preparing a surface

Two preprocessing steps, both optional and both pipeable:

```r
dem |> rvt_fill() |> rvt_smooth() |> rvt_curvature() |> rvt_plot()
```

`rvt_fill()` interpolates across NoData holes. Worth doing first, because a
hole is not just a missing answer — the derivative metrics lose a ring of cells
around every one, and the horizon metrics treat a hole as *transparent*, so a
gappy DSM reads brighter than the real surface.

`rvt_smooth()` removes noise while leaving breaks of slope sharp, using Sun et
al.'s feature-preserving method (as in WhiteboxTools' `FeaturePreservingSmoothing`).
That distinction is the whole point: on a noisy scarp a Gaussian of comparable
strength removes slightly more noise but smears the step over 6 cells, while
this keeps it exactly 1 cell wide. Useful before `rvt_curvature()`,
`rvt_log()` and `rvt_geomorphons()`, which all amplify high-frequency noise by
construction.

Fill before you smooth — smoothing cannot sensibly average across a hole.

### Scale of analysis

Most of these metrics answer a different question at a different resolution —
metre-scale detail on a 0.25 m DTM, landform structure on a 2 m version of the
same ground. `rvt_resample()` coarsens a raster (it refuses to go finer; detail
that was never measured can't be invented) and pipes straight on:

```r
dem |> rvt_resample(2) |> rvt_curvature() |> rvt_plot()
```

Use `method = "mode"` for categorical rasters such as `rvt_geomorphons()`
output, so classes stay classes.

Coarsening isn't the only scale knob, and often not the best one — for
`rvt_slrm()`, `rvt_msrm()`, `rvt_dev()` and `rvt_mstp()` a wider radius costs
nothing and keeps the fine detail in the input. Coarsening earns its place
where the metric has no radius at all (`rvt_slope()`, `rvt_aspect()`,
`rvt_curvature()`, `rvt_log()`), or to make a horizon-search metric affordable
— on a 4000×4000 0.25 m tile, sky-view factor over a 10 m reach takes 15.8 s,
against 0.6 s on the 1 m version of the same ground. `?rvt_resample` has the
full comparison.

### Looking a long way out

The horizon-search metrics — `rvt_svf()`, `rvt_asvf()`, both `rvt_openness()`
functions, `rvt_daylight()`, `rvt_shadow()` and `rvt_sky_illumination()` — scan
out to `reach` metres. Doing all of that at full resolution would make cost
grow in proportion to the distance, which is why "far enough to catch the ridge
across the valley" would be unaffordable on fine data.

So the full-resolution scan stops after 100 cells and everything beyond that is
read from coarser and coarser copies of the same DEM. Because 100 *cells* is
100 m on a 1 m DEM but only 25 m on a 0.25 m one, that is where the switch
happens: a `reach` inside it is a single, exact, full-resolution scan, and
anything longer adds levels — on a 0.25 m DEM, `reach = 1600` uses 0.25 m data
out to 25 m, then 1 m to 100 m, 4 m to 400 m and 16 m to 1600 m.

```r
dem |> rvt_svf(reach = 400)     # 400 m out, without a 400 px search
```

Cost then grows with the *logarithm* of the distance. Accuracy barely moves —
against a true full-resolution search of the same reach, sky-view factor differs
by 2×10⁻⁵ and openness by 0.0025°, because a feature's effect on the horizon is
`atan(height / distance)` and distant detail is worth very little. `?rvt_reach`
has the details, including the two knobs (`pyramid_px`, `pyramid_factor`) and
why the coarse levels use maximum rather than average resampling.

One caveat if your surface is a canopy or vegetation model rather than bare
earth: maximum resampling takes the tallest crown in each coarse block and
treats it as solid, so a lone tree or a mast blocks far more sky than it
really does. Pass `pyramid_method = "q3"` there — measured on a synthetic
canopy it halves the error, and the two are biased in opposite directions.

It matters when distant terrain genuinely shades or encloses the site. On flat
ground the far field contributes almost nothing, so a `reach` short enough to
stay within the full-resolution scan — and therefore exact — is all you need.

### Large rasters

Processing is tiled automatically; only one tile is in memory at a time, so input
size is limited by disk rather than RAM. Override the tile edge if you want:

```r
rvt_svf("huge_dem.tif", "svf.tif", tile_size = 4096, progress = TRUE)
```

Outputs are written as **Cloud-Optimized GeoTIFFs** — internally tiled, with
overviews built in. That means reading a result back doesn't have to mean
reading all of it: a viewer (including `rvt_plot()`) can ask for a fast
low-resolution overview instead of decoding the full raster, or a full-resolution
window over just part of a huge one — either way, GDAL only touches the bytes
that window actually needs.

### Threads

Defaults to all cores. `rvt_threads(8)` or `options(rvtr.threads = 8)` to change.

### Plotting

`rvt_plot()` gives a quick look at a result without a separate GIS viewer.
Downsampling happens on read via GDAL, so even a huge raster renders fast:

```r
dem |> rvt_vat("vat.tif") |> rvt_plot()
```

## Correctness

This isn't aiming for bit-for-bit parity with rvt-py — where rvt-py has a bug,
this package does the mathematically correct thing instead (see below).
`testthat::test_local()` checks, using the `dtm1.tif` tile bundled in
`inst/extdata/`:

- every metric above closely tracks rvt-py's output for the same DEM and
  default parameters
- a mosaic of tiles gives bit-identical results to the same area as one raster
- tile size has no effect on output (checked separately for the horizon-search
  kernel and the integral-image kernel behind MSRM, since they're different
  code paths)

All thirteen of RVT's visualization products are implemented.

Sky illumination searches the horizon at full resolution rather than through
RVT's multi-level DEM pyramid. The pyramid exists to make the default 100 px
radius affordable in numpy; the C++ kernel here is fast enough not to need the
approximation, at the cost of taking about a minute on a 4000 × 4000 raster —
drop `reach` if that matters more than exactness. With rvt-py forced to a
single resolution too, the two agree to a correlation of 0.997.

That gives confidence the algorithms are implemented correctly without locking
behaviour to another implementation's quirks. The rvt-py comparison outputs
(`inst/extdata/rvt_*.tif`) are **not** committed to git — they're a local
development aid for spot-checking against rvt-py, not a shipped fixture, so
those specific checks skip themselves when the files aren't present. See
`inst/extdata/README.md` for how to regenerate them.

Several rvt-py quirks are deliberately **not** reproduced by default:

- rvt-py's anisotropic sky-view factor brightens the sky at the **mirror** of
  the azimuth you ask for: its own default of 315° (north-west, matching its
  hillshade) actually lights from 45°, north-east. Its direction weights are
  computed from an angle that runs opposite to the direction the search
  actually steps in, and north/south happen to be fixed points of that mirror,
  so it only shows up east-west. `rvt_asvf(main_direction = 315)` here really
  does brighten the north-west, consistent with `rvt_hillshade()`. Verified on
  a synthetic scene with a wall to the east: the two implementations return
  the same numbers with east and west swapped. Sky illumination weights each
  direction against the terrain's aspect through the same mirrored angle, and
  is corrected here in the same way.
- **rvt-py's sky illumination cannot run on a DEM with NoData at all.** It
  pads its pyramid with `np.min(dem)`, which is `NaN` as soon as one cell is
  NoData, and that then propagates through every level — the output is
  entirely NaN. `rvt_sky_illumination()` handles NoData like every other
  metric here.
- rvt-py scales its overcast sky illumination by the brightest pixel in the
  raster. That cannot survive tiling — the same ground would get a different
  value depending on which tile it landed in, breaking the guarantee that
  results don't depend on how a raster is split. This package divides by the
  value flat ground produces instead, so flat reads 1 as it does for
  sky-view factor and local dominance.

- rvt-py treats "no horizon found in this direction" (only reachable beside
  NoData) as a slope of -1000 rather than -∞, off by a fraction of a degree in
  positive openness near NoData holes. This package uses the exact value.
- rvt-py's `blend_overlay()` mutates its `background` argument in place, so the
  VAT openness layer's declared 50% opacity ends up having no effect on its
  output. This package applies the opacity as configured. Pass
  `rvt_vat(rvt_compat = TRUE)` to reproduce rvt-py's output instead, e.g. to
  check results against it directly.

RVT's other padding choice is kept because it's the right one regardless of RVT:
reflect-pad for the horizon search (extrapolates terrain more plausibly at the
raster edge), edge-replicate for the slope/aspect derivatives (standard for
finite differences).
