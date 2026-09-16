# Which Mapterhorn sources cover a place

Lists the elevation sources Mapterhorn holds for an area, finest first,
without fetching any elevation data - to see what resolution
[`rvt_data_mapterhorn()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_mapterhorn.md)
will return, and whom to credit.

## Usage

``` r
rvt_data_mapterhorn_sources(x, size = 1000)
```

## Arguments

- x:

  location: a WGS84 `c(lon, lat)` point or `c(xmin, ymin, xmax, ymax)`
  extent, an sf point, geometry or bbox, or the path to a raster in a
  projected coordinate system

- size:

  side of the square returned for a point, in metres (default 1000)

## Value

a data frame with columns `source`, `name`, `producer`, `license`,
`resolution` (metres) and `website`

## See also

[`rvt_data_mapterhorn()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_mapterhorn.md)

## Examples

``` r
# \donttest{
rvt_data_mapterhorn_sources(c(9.9464, 51.6317))
#>            source                              name
#> 1 deniedersachsen Digitales Geländemodell 1m (DGM1)
#> 2           glo30                 COPERNICUS GLO-30
#>                                                                                                                                                                                    producer
#> 1                                                                                                                  Landesamtes für Geoinformation und Landesvermessung Niedersachsen (LGLN)
#> 2 DLR e.V. 2010-2014, Airbus Defence and Space GmbH 2014-2018, provided under COPERNICUS by the European Union and ESA. All rights reserved. Accessed via OpenTopography and AWS Open Data.
#>                                  license resolution
#> 1                              CC BY 4.0          1
#> 2 COPERNICUS full, free and open license         30
#>                                                                                                                          website
#> 1      https://www.lgln.niedersachsen.de/startseite/geodaten_karten/3d_geobasisdaten/dgm/digitale-gelandemodelle-dgm-143150.html
#> 2 https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM
# }
```
