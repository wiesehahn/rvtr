# Geomorphon landform classes

The ten classes
[`rvt_geomorphons()`](https://wiesehahn.github.io/rvtr/reference/rvt_geomorphons.md)
assigns, in code order: a named integer vector mapping class name to the
value stored in the raster.

## Usage

``` r
rvt_geomorphon_classes
```

## Format

Named integer vector of length 10.

## Examples

``` r
rvt_geomorphon_classes
#>      flat      peak     ridge  shoulder      spur     slope    hollow footslope 
#>         1         2         3         4         5         6         7         8 
#>    valley       pit 
#>         9        10 
names(rvt_geomorphon_classes)[3]
#> [1] "ridge"
```
