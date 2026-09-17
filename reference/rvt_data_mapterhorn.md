# Terrain for anywhere in the world, from Mapterhorn

A digital terrain model for any place on Earth, from
[Mapterhorn](https://mapterhorn.com), which combines open elevation data
into one global dataset: 30 m Copernicus GLO-30 everywhere, and national
LiDAR terrain models where they have been published - 1 m across
Germany, Austria and Norway, 0.5 m in Switzerland, 0.25 m in parts of
Baden-Württemberg, and many more. Give a point to get a square around
it, or an extent to get exactly that area.

## Usage

``` r
rvt_data_mapterhorn(
  x,
  res = NULL,
  size = 1000,
  crs = NULL,
  download = FALSE,
  max_mb = 500,
  refresh = FALSE,
  threads = rvt_threads()
)
```

## Arguments

- x:

  location: a WGS84 `c(lon, lat)` point or `c(xmin, ymin, xmax, ymax)`
  extent, an sf point, geometry or bbox, or the path to a raster in a
  projected coordinate system

- res:

  output resolution in metres, or `NULL` (default) for the finest source
  available

- size:

  side of the square returned for a point, in metres (default 1000)

- crs:

  output coordinate system: an EPSG code, `"EPSG:..."` string or WKT, in
  metres. `NULL` (default) uses the local UTM zone.

- download:

  `FALSE` (default) returns a virtual raster reading the tiles on
  demand; `TRUE` stores a local copy in the user cache

- max_mb:

  stop before reading if the estimated download exceeds this many
  megabytes (default 500)

- refresh:

  with `download = TRUE`, fetch again even if cached

- threads:

  C++ threads (see
  [`rvt_threads()`](https://wiesehahn.github.io/rvtr/reference/rvt_threads.md))

## Value

path to the raster (a VRT, or a Cloud-Optimized GeoTIFF with
`download = TRUE`), so it pipes into any function

## Details

Nothing is downloaded by default. The result is a small virtual raster
that reads Mapterhorn's tiles only when something reads it, and only the
tiles that read needs. Pass it to any function in this package.
`download = TRUE` stores a local copy instead, for offline work or
repeated sessions.

## Resolution

`res = NULL` uses the finest source covering the area; see what that is
with
[`rvt_data_mapterhorn_sources()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_mapterhorn_sources.md).
Asking for a finer `res` than that is an error, because it could only
interpolate. Where the area straddles sources of different resolution -
the edge of a national LiDAR survey, say - the coarser source fills in
wherever the finer one has no data, so the result has no holes.

Heights are stored to 1/256 m (about 4 mm), which is below the accuracy
of any source.

## Coordinate system

Mapterhorn's tiles are in Web Mercator, whose map units are stretched by
1/cos(latitude): 1.6 times at 51 degrees north. Since every distance in
this package is in map units, the data is always reprojected to a
coordinate system in metres - by default the UTM zone of the area's
centre (standard zones, e.g. EPSG:32632), or whatever `crs` you pass.
Pass the national system to line data up with other sources, e.g.
`crs = 25832` for Germany.

## Location

- **A point** returns a `size` x `size` m square centred on it. Plain
  numbers `c(lon, lat)` are WGS84; an sf point uses its own coordinate
  system.

- **An extent** returns exactly that area, snapped outward to whole
  cells: an
  [`sf::st_bbox()`](https://r-spatial.github.io/sf/reference/st_bbox.html),
  any sf geometry, or plain WGS84 `c(xmin, ymin, xmax, ymax)`.

- **A raster path** returns data on exactly that raster's grid and
  coordinate system, which must be projected.

## Licence

Mapterhorn's data comes from many producers under their own open
licences (CC BY 4.0, Datenlizenz Deutschland, and others). The sources
used are listed in a message the first time each is fetched in a
session, and by
[`rvt_data_mapterhorn_sources()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_mapterhorn_sources.md).
Credit Mapterhorn and each source; see
<https://mapterhorn.com/attribution> for the full terms.

## See also

[`rvt_data_mapterhorn_sources()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_mapterhorn_sources.md),
and
[`rvt_data_lgln()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln.md)
for Lower Saxony's terrain, surface model and orthophotos at their
native grid

## Examples

``` r
# \donttest{
# the Matterhorn, 3 km across at Switzerland's 0.5 m, shown at 2 m
rvt_data_mapterhorn(c(7.6586, 45.9763), res = 2, size = 3000) |>
  rvt_hillshade() |>
  rvt_plot()
#> Mapterhorn terrain, https://mapterhorn.com/attribution. Sources here:
#>   swissALTI3D - Federal Office of Topography swisstopo, Open Government Data
#>   Modello Digitale del Terreno (DTM) - Dati estratti dal Modello Digitale del Terreno (DTM) della Regione Autonoma Valle d'Aosta, CC BY 4.0
#>   TINITALY, a digital elevation model of Italy with a 10 meters cell size (Version 1.1) - Istituto Nazionale di Geofisica e Vulcanologia (INGV), Creative Commons Namensnennung 4.0 International (CC BY 4.0)
#>   COPERNICUS GLO-30 - DLR e.V. 2010-2014, Airbus Defence and Space GmbH 2014-2018, provided under COPERNICUS by the European Union and ESA. All rights reserved. Accessed via OpenTopography and AWS Open Data., COPERNICUS full, free and open license
#> Warning: GDAL Error 5: VRTDerivedRasterBand::IRasterIO:Derived band pixel function 'expression' not registered.
#> Warning: GDAL Error 1: /tmp/RtmpXpCjHd/file218c221a3081.vrt, band 1: IReadBlock failed at X offset 0, Y offset 0: VRTDerivedRasterBand::IRasterIO:Derived band pixel function 'expression' not registered.
#> Error: read raster failed

rvt_data_mapterhorn_sources(c(7.6586, 45.9763))
#>        source
#> 1 swissalti3d
#> 2     itaosta
#> 3    tinitaly
#> 4       glo30
#>                                                                                    name
#> 1                                                                           swissALTI3D
#> 2                                                    Modello Digitale del Terreno (DTM)
#> 3 TINITALY, a digital elevation model of Italy with a 10 meters cell size (Version 1.1)
#> 4                                                                     COPERNICUS GLO-30
#>                                                                                                                                                                                    producer
#> 1                                                                                                                                                    Federal Office of Topography swisstopo
#> 2                                                                                                 Dati estratti dal Modello Digitale del Terreno (DTM) della Regione Autonoma Valle d'Aosta
#> 3                                                                                                                                     Istituto Nazionale di Geofisica e Vulcanologia (INGV)
#> 4 DLR e.V. 2010-2014, Airbus Defence and Space GmbH 2014-2018, provided under COPERNICUS by the European Union and ESA. All rights reserved. Accessed via OpenTopography and AWS Open Data.
#>                                                        license resolution
#> 1                                         Open Government Data        0.5
#> 2                                                    CC BY 4.0        2.0
#> 3 Creative Commons Namensnennung 4.0 International (CC BY 4.0)       10.0
#> 4                       COPERNICUS full, free and open license       30.0
#>                                                                                                                          website
#> 1                                                                     https://www.swisstopo.admin.ch/en/height-model-swissalti3d
#> 2                                                                                https://geoportale.regione.vda.it/download/dtm/
#> 3                                                                              https://tinitaly.pi.ingv.it/Download_Area1_1.html
#> 4 https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM
# }
```
