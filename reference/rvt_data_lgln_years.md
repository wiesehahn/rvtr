# Years of Lower Saxony data available for a place

Lists the flight years LGLN holds for a location, without fetching any
data: one row per product and year, with the flight dates, the number of
tiles involved, and whether that year covers the whole area.

## Usage

``` r
rvt_data_lgln_years(x, product = c("dtm", "dsm", "rgb"))
```

## Arguments

- x:

  location: a WGS84 `c(lon, lat)` point or `c(xmin, ymin, xmax, ymax)`
  extent, an sf point, geometry or bbox, or the path to an EPSG:25832
  raster

- product:

  one or more of `"dtm"`, `"dsm"`, `"rgb"` (default all)

## Value

a data frame with columns `product`, `year`, `first_date`, `last_date`,
`tiles` and `covers`

## See also

[`rvt_data_lgln()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln.md)

## Examples

``` r
# \donttest{
rvt_data_lgln_years(c(9.9464, 51.6317))
#>   product year first_date  last_date tiles covers
#> 1     dtm 2016 2016-04-01 2016-04-01     1   TRUE
#> 2     dsm 2016 2016-04-01 2016-04-01     1   TRUE
#> 3     rgb 2013 2013-06-18 2013-06-18     1   TRUE
#> 4     rgb 2016 2016-03-17 2016-03-17     1   TRUE
#> 5     rgb 2019 2019-04-01 2019-04-01     1   TRUE
#> 6     rgb 2022 2022-03-13 2022-03-13     1   TRUE
#> 7     rgb 2025 2025-03-08 2025-03-08     1   TRUE
# }
```
