# Curvature

How sharply the ground bends, and in which sense: **positive is convex**
(the crest of a bank, the lip of a terrace), **negative concave** (the
floor of a ditch, the foot of a slope), zero where the surface is
locally planar - whether that plane is flat or steeply tilted.

## Usage

``` r
rvt_curvature(
  dem,
  out_path = fs::file_temp(ext = "tif"),
  type = c("profile", "plan", "tangential", "mean", "minimal", "maximal", "gaussian",
    "total", "unsphericity", "casorati", "shape_index", "difference", "twisting",
    "rotor"),
  radius = NULL,
  tile_size = NULL,
  threads = rvt_threads(),
  overwrite = FALSE,
  progress = FALSE
)
```

## Arguments

- dem:

  path to a DEM, or a vector of paths forming a mosaic (see
  [`rvt_mosaic()`](https://wiesehahn.github.io/rvtr/reference/rvt_mosaic.md));
  mosaics are read across file boundaries, so tiles do not produce
  seams. A terra `SpatRaster` works too: one read from a file is used as
  it lies, while one terra computed in memory is written to a temporary
  GeoTIFF first.

- out_path:

  output GeoTIFF path (default: a temp file). A folder that does not
  exist yet is created.

- type:

  one of `"profile"` (default), `"plan"`, `"tangential"`, `"mean"`,
  `"minimal"`, `"maximal"`, `"gaussian"`, `"total"`, `"unsphericity"`,
  `"casorati"`, `"shape_index"`, `"difference"`, `"twisting"`, `"rotor"`

- radius:

  half-width of the window a quadratic surface is fitted over, in **map
  units**; omit it to measure across single pixels instead

- tile_size:

  tile edge in pixels; NULL picks a size targeting roughly 256 MB per
  working matrix

- threads:

  C++ threads (see
  [`rvt_threads()`](https://wiesehahn.github.io/rvtr/reference/rvt_threads.md))

- overwrite:

  recompute even if `out_path` exists

- progress:

  report per-tile progress

## Value

`out_path`, invisibly

## Details

Curvature is the natural detector for *breaks of slope*, which is often
what actually defines an earthwork: a bank is a convex line beside a
concave one. It is also the sharpest of the metrics here, in both
senses - it responds to the finest detail in the DEM, and it responds
just as eagerly to noise. Expect to need a smoothed DEM, or to prefer
[`rvt_log()`](https://wiesehahn.github.io/rvtr/reference/rvt_log.md)
which folds the smoothing in.

Units are per map unit (1/m on a metric DEM), and values are small: a
subtle bank might reach 0.05. There is no natural range, so displaying
it means picking a symmetric stretch around zero - something like -0.1
to 0.1

- and a diverging palette so convex and concave read as opposite.
  `shape_index` is the exception: it is dimensionless and already
  bounded to -1 to 1, so it displays as it comes.

## How it works

A quadratic surface (Zevenbergen and Thorne 1987) is fitted to the eight
neighbours of each cell, and the named curvatures are combinations of
its second derivatives:

- `profile`:

  **Does water speed up or slow down here?** Positive where a slope
  steepens towards the brow of a bank, negative where it eases off at
  the foot. Measured straight down the slope. This is the one that picks
  out breaks of slope, and usually what people mean by "curvature"
  unqualified.

- `plan`:

  **Does water gather together or spread apart here?** Negative in
  hollows and gullies where flow converges, positive on spurs and noses
  where it fans out. Measured across the slope, along the contour.

- `tangential`:

  The same gathering-or-spreading question as `plan`, but better behaved
  on gentle ground, where contours wander and plan curvature turns
  erratic. Prefer it on anything flat.

- `mean`:

  **Does the ground bulge outwards or dish inwards here?** Averaged over
  every direction at once. The general-purpose convexity measure, and
  the only one of these that stays meaningful on level ground.

- `minimal`,`maximal`:

  **Which way does the ground bend most, and which way least?** A ridge
  bends hard across its crest and barely at all along it, so `maximal`
  is large and positive while `minimal` sits near zero; a valley is the
  mirror of that. A dome has both positive, a pit both negative.

- `gaussian`:

  **Is this a dome or bowl, or a saddle?** Positive where the ground
  bends the same way in every direction, negative where it bends up one
  way and down the other - passes between two summits, or the junction
  of two valleys.

- `total`:

  **How much is happening here at all?** How strongly the ground bends
  in any direction, ignoring whether it bends up or down. A magnitude,
  so never negative - useful for finding where something is going on
  when you do not care what.

- `unsphericity`:

  **Is the bending even or lopsided?** Near zero where the ground curves
  the same way in every direction - a dome, a bowl, a plane - and large
  on ridges, valleys and saddles that bend mostly one way. Never
  negative.

- `casorati`:

  **How sharp is the bending, whatever shape it makes?** Zero only where
  the ground is genuinely flat. Like `total` in spirit but a properly
  normalised curvature, so it stays comparable between gentle and steep
  ground.

- `shape_index`:

  **What shape is the ground, regardless of how hard it bends?** A pit
  reads -1, a valley -0.5, a saddle 0, a ridge +0.5 and a peak +1. It
  separates form from intensity, so a faint ridge and a sharp one read
  alike, and being bounded it needs no stretch to display.

- `difference`:

  **Is this place more about water speeding up, or more about water
  gathering?** Mean curvature averages those two and can read zero where
  a strong acceleration cancels a strong convergence; this is the
  distinction that goes missing there.

- `twisting`:

  **Is the hillside twisting?** How fast the steepness changes as you
  walk *sideways* across a slope rather than up or down it. Zero on any
  slope of even steepness, however steep. Positive where the ground is
  steeper to your right as you look downhill. It is exactly the rate of
  change of slope angle along the contour, in radians per map unit, so
  the number means something directly.

- `rotor`:

  **Does the downhill path bend as it descends?** Positive where a ball
  rolling down would veer right, negative left, zero where it runs
  straight. Picks out corkscrewing spurs and gullies that turn as they
  fall - a thing none of the others above can see, because they all ask
  about bending *across* or *along* the slope, never about the slope
  direction itself changing. It is the projected form of `twisting` and
  shares `plan`'s weakness for the same reason: on near-level ground the
  downhill direction is barely defined, so values explode (hundreds,
  against `twisting`'s fraction of one). Use `twisting` unless you
  specifically want the plan-view bend.

`unsphericity`, `casorati`, `shape_index` and `difference` come from
Shary's system, of which `mean`, `unsphericity` and `difference` are the
three independent components - every other curvature here is a
combination of those. `twisting` completes the basic trio of `profile`,
`plan` and twisting, and `rotor` is its projected form.

Except for `total`, these are the proper differential-geometry
quantities, carrying the `(1 + slope^2)` denominators that turn second
derivatives into actual curvature. Dropping those - as some
implementations do - agrees only where the ground is level and
overstates curvature increasingly with slope, by more than a factor of
two by 17 degrees. `total` is the exception on purpose: it follows
Evans' un-normalised definition, because that is what the name means in
the geomorphometry literature and what other software returns for it.

Checked against closed-form values on a tilted paraboloid and against
WhiteboxTools: `mean`, `gaussian`, `minimal`, `maximal` and `profile`
agree with both to four significant figures, and `total` matches
WhiteboxTools.

Give a `radius` and the scale changes: a quadratic surface is fitted by
least squares over the whole window (Evans 1980; Wood 1996) and the same
curvatures are read off its coefficients. This is the standard way to
pick a scale, and the only one available - the Zevenbergen and Thorne
fit passes exactly through nine points and has no wider form. The two
are not identical at the smallest radius: the fit averages the bending
over every row of the window, where the default takes the middle row
alone. Leaving `radius` out therefore keeps the values other software
returns, and setting it opts into the multi-scale form.

Profile and plan curvature need a slope direction to be defined at all;
on genuinely flat cells they are reported as 0 rather than as a division
by nearly nothing. NoData inside the window contributes the centre
cell's own elevation, so a hole never eats a ring of cells around
itself - but at a large `radius` that biases curvature towards zero near
big holes, so fill them first with
[`rvt_fill()`](https://wiesehahn.github.io/rvtr/reference/rvt_fill.md).

## Tuning

- `type` is the only real choice, and `"profile"` is the usual starting
  point for finding earthworks.

- `radius` sets the scale, in **map units**. Omitted, curvature is
  measured across single pixels, which on a noisy DEM shows. A few
  metres steadies it and picks out banks and ditch lips; tens of metres
  describes the shape of the hillslope instead. Cost grows with the
  square of the radius, so keep it to what the features need.

- [`rvt_log()`](https://wiesehahn.github.io/rvtr/reference/rvt_log.md)
  is the alternative when you want an edge detector rather than a
  curvature: its `sigma` sets a scale too, but by smoothing first.

## References

Zevenbergen, L. W. and Thorne, C. R. (1987) Quantitative analysis of
land surface topography. *Earth Surface Processes and Landforms* 12,
47-56.

Evans, I. S. (1980) An integrated system of terrain analysis and slope
mapping. *Zeitschrift für Geomorphologie* Supplementband 36, 274-295.

Wood, J. (1996) *The Geomorphological Characterisation of Digital
Elevation Models*. PhD thesis, University of Leicester.

Shary, P. A. (1995) Land surface in gravity points classification by a
complete system of curvatures. *Mathematical Geology* 27, 373-390.
[doi:10.1007/BF02084608](https://doi.org/10.1007/BF02084608)

Koenderink, J. J. and van Doorn, A. J. (1992) Surface shape and
curvature scales. *Image and Vision Computing* 10, 557-564.
[doi:10.1016/0262-8856(92)90076-F](https://doi.org/10.1016/0262-8856%2892%2990076-F)

Minár, J., Evans, I. S. and Jenčo, M. (2020) A comprehensive system of
definitions of land surface (topographic) curvatures. *Earth-Science
Reviews* 211, 103414.
[doi:10.1016/j.earscirev.2020.103414](https://doi.org/10.1016/j.earscirev.2020.103414)

## See also

[`rvt_log()`](https://wiesehahn.github.io/rvtr/reference/rvt_log.md) for
an edge detector with a scale parameter.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")
dem |> rvt_curvature()
dem |> rvt_curvature(type = "plan")
dem |> rvt_curvature(radius = 5)
```
