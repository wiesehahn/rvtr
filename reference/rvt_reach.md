# Reaching further with a multi-resolution search

The horizon-search metrics -
[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md),
[`rvt_asvf()`](https://wiesehahn.github.io/rvtr/reference/rvt_asvf.md),
[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md),
[`rvt_openness_negative()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness_negative.md),
[`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md),
[`rvt_shadow()`](https://wiesehahn.github.io/rvtr/reference/rvt_shadow.md)
and
[`rvt_sky_illumination()`](https://wiesehahn.github.io/rvtr/reference/rvt_sky_illumination.md) -
scan outwards a fixed distance, `reach`. Scanning all of it at full
resolution makes the cost grow in proportion to how far you look, which
is why a reach long enough to catch the ridge across the valley would be
unaffordable on fine data.

So they don't scan all of it at full resolution. The full-resolution
scan stops after `pyramid_px` cells - 100 by default - and everything
beyond that is read from coarser and coarser copies of the same DEM,
each `pyramid_factor` times coarser than the last (4 by default) and
each again `pyramid_px` cells wide. Cost then grows with the *logarithm*
of the reach instead of in proportion to it.

## What that means for a given DEM

That first 100 is a count of **cells**, so how far it reaches depends on
the resolution: 100 m on a 1 m DEM, 25 m on a 0.25 m one, 200 m on a 2 m
one. Ask for a `reach` within that and the search is a single
full-resolution scan - exactly what it has always been, and exact. Ask
for more and levels are added behind the scenes:

    0.25 m DEM, reach   20 m  ->  0.25 m data throughout
    0.25 m DEM, reach  100 m  ->  0.25 m out to 25 m, then 1 m to 100 m
    0.25 m DEM, reach 1600 m  ->  0.25 m out to 25 m, 1 m to 100 m,
                                  4 m to 400 m, 16 m to 1600 m
       1 m DEM, reach  400 m  ->  1 m out to 100 m, then 4 m to 400 m
       2 m DEM, reach  100 m  ->  2 m data throughout

Nothing else about the call changes, and the coarse copies are built,
used and discarded for you.

## Why it works

A feature's effect on the horizon is an angle,
`atan(height / distance)`, so detail matters in proportion to how close
it is. A 1 m bump 5 m away is 11 degrees of sky; the same bump 400 m
away is 0.14 degrees. Sampling the far field on a coarse grid therefore
costs almost nothing in accuracy, while saving almost all of the work.

Level 0 is the DEM itself, searched to `pyramid_px` pixels. Each further
level is `pyramid_factor` times coarser and covers from where the last
one stopped out to `pyramid_factor` times further, again in `pyramid_px`
pixels. Each cell's horizon is the steepest angle found across all
levels, always measured from that cell's own full-resolution elevation.

The coarse levels are built with **maximum** resampling, not averaging.
Averaging erases anything narrower than a coarse cell - exactly the
isolated crags, walls and tree lines that make up a distant skyline. A 3
m wall 20 m high at 300 m, coarsened 16 times, keeps 3.47 of its 3.53
degrees of horizon under `"max"`, against 0.65 under cubic and 0.00
under nearest.

## When max is the wrong choice

On a **rough** surface that is all peaks - a canopy or vegetation
model - max errs the other way: it takes the tallest crown in each block
and treats it as solid ground, so a lone tree or a mast is widened to
the whole coarse cell and blocks far more sky than it really does.
Measured on a synthetic canopy, mean sky-view error against a
full-resolution search is 0.0066 for `"max"` against **0.0030 for
`"q3"`**, and the two are biased opposite ways: max over-states
obstruction, the third quartile slightly under-states it.

So `pyramid_method = "q3"` on canopy surfaces, and the default `"max"`
on bare earth, where q3 would erase a real wall entirely. A feature has
to occupy about a quarter of a coarse cell to survive q3, which is
roughly the point at which it starts blocking a meaningful share of that
cell's sky.

Negative openness searches inverted terrain, so its levels use the
mirror-image statistic automatically - min for max, q1 for q3 - since
the horizon it wants is `max(-z)`, which is `-min(z)`.

## What it costs in accuracy

Measured against a full-resolution search of the same reach, on terrain
with 300 m of relief where the 100-400 m band is worth 0.74 degrees of
openness, the defaults are in error by 0.068 degrees - about 9% of the
far-field signal they are approximating, and far less than that signal
itself. On gentle terrain the error is nearer 0.003 degrees.

`pyramid_px` is the accuracy knob, and it matters much more than
`pyramid_factor`: doubling it roughly halves the error, because it is
what decides how far *full resolution* extends. Raising `pyramid_factor`
is the cheaper, blunter lever.

## When to use it

Use `reach` when distant terrain genuinely shades or encloses the site -
mountains, deep valleys, quarry faces - or when you want a long reach at
fine resolution and cannot afford it otherwise. On flat or gently
rolling ground the far field contributes very little (0.0004 of sky-view
factor on this package's sample tile), so a `reach` short enough to stay
within the full-resolution scan - and therefore exact - is all you need.

`reach` is in **map units** (metres, normally), as every distance in
this package is, so it means the same thing whatever the resolution of
the DEM. To force an exact full-resolution search at any distance, raise
`pyramid_px` until one level covers the whole reach.

## See also

[`rvt_svf()`](https://wiesehahn.github.io/rvtr/reference/rvt_svf.md),
[`rvt_openness()`](https://wiesehahn.github.io/rvtr/reference/rvt_openness.md),
[`rvt_daylight()`](https://wiesehahn.github.io/rvtr/reference/rvt_daylight.md),
and
[`rvt_resample()`](https://wiesehahn.github.io/rvtr/reference/rvt_resample.md),
which builds the coarse levels.

## Examples

``` r
dem <- system.file("extdata", "dtm1.tif", package = "rvtr")

# look 400 m out from a 1 m DEM without paying for a 400 px search
dem |> rvt_svf(reach = 400)
```
