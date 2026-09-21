# Years of Lower Saxony data available for a place

Lists the flight years LGLN holds for a location, without fetching any
data: one row per product and year, with the flight dates, the number of
tiles involved, and how much of the area asked for that year's tiles
fill.

## Usage

``` r
rvt_data_lgln_years(x, product = c("dtm", "dsm", "bdom", "rgb"))
```

## Arguments

- x:

  location: a WGS84 `c(lon, lat)` point or `c(xmin, ymin, xmax, ymax)`
  extent, an sf point, geometry or bbox, or the path to an EPSG:25832
  raster

- product:

  one or more of `"dtm"`, `"dsm"`, `"bdom"`, `"rgb"` (default all)

## Value

a data frame, one row per product and year:

- `product`:

  `"dtm"`, `"dsm"`, `"bdom"` or `"rgb"`

- `year`:

  the flight year

- `first_date`, `last_date`:

  earliest and latest flight date among that year's tiles, as
  `"YYYY-MM-DD"`. They differ when neighbouring tiles were flown on
  different days of the same campaign

- `tiles`:

  how many distinct tiles that year contributes to the area

- `coverage`:

  the fraction of the area those tiles fill, 0 to 1 - see above. Below
  1, naming that year in
  [`rvt_data_lgln()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln.md)
  returns that share and leaves the rest NoData

- `covers`:

  `coverage == 1`, as a convenience: the years that can be named on
  their own and still fill the whole area

## What `coverage` and `covers` mean

The catalogue is a patchwork of 1 km and 2 km tiles, and any one year is
flown over only part of Lower Saxony. `coverage` is the exact fraction
of your area that year's tiles fill, and `covers` is the same thing as a
yes-or-no (`coverage` of 1). Anything below 1 means an area lying across
the boundary between two flight campaigns, where that year's tiles stop
partway across it.

It is what to read before naming a `year` in
[`rvt_data_lgln()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln.md),
because a named year is used on its own: a year at 0.43 gives you that
43% and NoData for the rest (or, under `partial = FALSE`, an error).
Leaving `year` unset instead gives each tile its own latest flight,
which covers as much ground as all these rows together - so a set of
rows that are each partial can still add up to a complete raster,
assembled from several years.

A point is always answered with the single 1 km tile containing it, so
`coverage` is 1 for every year listed. Only an extent, an sf bbox or a
raster grid can lie across a seam.

## See also

[`rvt_data_lgln()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln.md)

## Examples

``` r
# \donttest{
rvt_data_lgln_years(c(9.9464, 51.6317))
#>   product year first_date  last_date tiles coverage covers
#> 1     dtm 2016 2016-04-01 2016-04-01     1        1   TRUE
#> 2     dsm 2016 2016-04-01 2016-04-01     1        1   TRUE
#> 3    bdom 2022 2022-03-13 2022-03-13     1        1   TRUE
#> 4    bdom 2025 2025-03-08 2025-03-08     1        1   TRUE
#> 5     rgb 2013 2013-06-18 2013-06-18     1        1   TRUE
#> 6     rgb 2016 2016-03-17 2016-03-17     1        1   TRUE
#> 7     rgb 2019 2019-04-01 2019-04-01     1        1   TRUE
#> 8     rgb 2022 2022-03-13 2022-03-13     1        1   TRUE
#> 9     rgb 2025 2025-03-08 2025-03-08     1        1   TRUE
# }
```
