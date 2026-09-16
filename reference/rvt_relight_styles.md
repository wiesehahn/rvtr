# Lighting styles for rvt_relight()

Named sets of lighting settings for
[`rvt_relight()`](https://wiesehahn.github.io/rvtr/reference/rvt_relight.md).
`default` holds every setting; each other style lists only what it
changes from `default`. Any setting passed to
[`rvt_relight()`](https://wiesehahn.github.io/rvtr/reference/rvt_relight.md)
directly overrides the style.

## Usage

``` r
rvt_relight_styles
```

## Format

A named list of lists.

## Details

- `default`:

  Warm afternoon sun from the north-west with soft blue shadows and
  slightly lifted colour.

- `neutral`:

  The same light without colour: grey sun and sky, to see what the tints
  add.

- `no_cast_shadows`:

  Shading and occlusion only; nothing casts a shadow onto its
  surroundings.

- `soft_shadows`:

  Lighter, more open shadows.

- `vivid`:

  Stronger colour and contrast.

- `low_sun`, `high_sun`:

  A lower sun for long shadows, or a higher one for short shadows.

- `hard_shadows`:

  Crisp shadow edges instead of soft ones.

- `verdant`:

  Woodland and fields brought forward: deeper, more varied greens, with
  everything else left as it is.

- `golden_hour`:

  Warm, low-contrast evening light.

- `winter`:

  Cold, frosty light with longer shadows and dimmer highlights.

- `twilight`:

  A very low rose-peach sun under a lavender sky.

- `blue_hour`:

  The sun barely present; the scene lit mostly by a cool sky.

## See also

[`rvt_relight()`](https://wiesehahn.github.io/rvtr/reference/rvt_relight.md)

## Examples

``` r
names(rvt_relight_styles)
#>  [1] "default"         "neutral"         "no_cast_shadows" "soft_shadows"   
#>  [5] "vivid"           "low_sun"         "high_sun"        "hard_shadows"   
#>  [9] "verdant"         "golden_hour"     "winter"          "twilight"       
#> [13] "blue_hour"      
rvt_relight_styles$golden_hour
#> $sun_color
#> [1] "#FFE0B8"
#> 
#> $sky_color
#> [1] "#8599D6"
#> 
#> $saturation
#> [1] 1.06
#> 
```
