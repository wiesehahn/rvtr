# Lower Saxony elevation and orthophotos, for any place

Terrain model, surface model or orthophoto for any location in Lower
Saxony, from the state survey authority (LGLN), in any flight year it
holds. Give a point to get the 1 km tile it lies in, an extent to get
exactly that area, or a raster to get data on its grid.

## Usage

``` r
rvt_data_lgln(
  x,
  product = c("dtm", "dsm", "bdom", "rgb"),
  year = NULL,
  res = NULL,
  download = FALSE,
  partial = TRUE,
  max_mb = 500,
  refresh = FALSE,
  threads = rvt_threads()
)
```

## Arguments

- x:

  location: a WGS84 `c(lon, lat)` point or `c(xmin, ymin, xmax, ymax)`
  extent, an sf point, geometry or bbox, or the path to an EPSG:25832
  raster

- product:

  `"dtm"` (default), `"dsm"`, `"bdom"` or `"rgb"`

- year:

  one flight year, giving a single survey throughout, or `NULL`
  (default) to give each tile its own latest flight - see the Years
  section, and
  [`rvt_data_lgln_years()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln_years.md)
  for what exists

- res:

  output resolution in metres, or `NULL` (default) for native - or, when
  `x` is a raster, for that raster's own resolution

- download:

  `FALSE` (default) returns a virtual raster reading the remote files on
  demand; `TRUE` stores a local copy in the user cache

- partial:

  `TRUE` (default) returns data wherever the chosen tiles reach and
  leaves the rest NoData; `FALSE` errors instead of returning a hole

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
that points at LGLN's published Cloud-Optimized GeoTIFFs, and data moves
only when something reads it - and only the parts that read needs. Pass
it to any function in this package. `download = TRUE` stores a local
copy instead, for offline work or repeated sessions.

|  |  |  |  |
|----|----|----|----|
| product | what | native grid | years |
| `"dtm"` | DGM1 terrain model, the bare ground | 1 m, 1 km tiles | usually one |
| `"dsm"` | DOM1 surface model, from lidar: ground plus trees and buildings | 1 m, 1 km tiles | usually one |
| `"bdom"` | bDOM20 surface model, matched from the aerial photos | 0.2 m, 1 km tiles | 2021 onwards |
| `"rgb"` | DOP20 orthophoto | 0.2 m, 2 km tiles | several, e.g. 2013-2025 |

Everything comes in ETRS89 / UTM zone 32N (EPSG:25832), the data's own
coordinate system; locations are transformed to find the data, but the
data is never reprojected.
[`rvt_data_lgln_years()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln_years.md)
lists which years exist for a place without fetching any data.

## Two surface models, and why the flight date matters

`"dsm"` is measured by lidar; `"bdom"` is computed by matching features
between overlapping aerial photographs. bDOM is five times finer - 0.2 m
against 1 m - and renders roofs, walls, hedges and open ground more
crisply, which is what
[`rvt_relight()`](https://wiesehahn.github.io/rvtr/reference/rvt_relight.md)
and the surface metrics want in a town.

**Matching only works where the photographs show matchable texture, so
over woodland the result depends on when the photographs were taken.** A
canopy in leaf gives plenty of texture and is captured well. A bare
broadleaf canopy gives almost none, and the match then falls through to
whatever does match, usually the forest floor - so bDOM can report *no
trees where trees stand*. That failure does not announce itself: the
wooded ground comes back smooth and plausible, and a hillshade of it
shows forest tracks and gullies as though the canopy were not there.

So check the season before trusting bDOM under trees.
[`rvt_data_lgln_years()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln_years.md)
gives the flight dates of every year available, and `year` picks one: a
summer flight may describe the canopy well, while the early-spring
flights that aerial survey usually favours will not. Where canopy
matters regardless - canopy height, forestry, relighting a wooded
scene - `"dsm"` is measured rather than inferred, and does not depend on
the season.

bDOM is also the heaviest product here by a distance: uncompressed 0.2 m
Float32, about 100 MB per km² at native resolution against the
orthophoto's 4.6 MB. Pass `res` unless you need every pixel; the cost
guard stops the call before it starts if you forget.

## Location

- **A point** returns the 1 x 1 km tile containing it. Plain numbers
  `c(lon, lat)` are WGS84; an sf point uses its own coordinate system,
  so EPSG:25832 coordinates work as an sf point.

- **An extent** returns exactly that area, snapped outward to whole
  cells: an
  [`sf::st_bbox()`](https://r-spatial.github.io/sf/reference/st_bbox.html),
  any sf geometry, or plain WGS84 `c(xmin, ymin, xmax, ymax)`.

- **A raster path** returns data on exactly that raster's grid, which
  must be EPSG:25832.

## Resolution and download size

`res = NULL` keeps native resolution. A number resamples, cubic, onto a
grid of whole multiples of `res`, so products fetched for the same place
at the same `res` line up exactly and can go straight into
[`rvt_blend()`](https://wiesehahn.github.io/rvtr/reference/rvt_blend.md).

Coarser requests read from the files' overviews, which cuts orthophoto
downloads roughly tenfold at 1 m. **The terrain and surface files have
no overview coarser than 2 m**, though, so even a 100 m request still
reads about 1 MB per km² - all of Lower Saxony would be around 50 GB.
Before anything is read, the download is estimated and the call stops
above `max_mb`. For large areas at coarse resolution, use a coarse
elevation source instead.

## Years

Any one year is flown over only part of Lower Saxony, so a single year's
tiles may fill your area or stop partway across it.

**`year = NULL` (the default) gives every tile its own latest flight**,
so you get the newest picture of the ground *and* the widest coverage:
where 2024 reaches you get 2024, and where it does not you get whatever
year last covered that tile. The cost is that neighbouring tiles can be
surveys years apart - trees grown, buildings put up - which can show as
a seam. A message names the years whenever more than one is used.

**Naming a `year` means that year and no other**, one survey throughout.
Where it was not flown you get NoData, and a message names the share
covered;
[`rvt_data_lgln_years()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln_years.md)
reports the same share up front, in its `coverage` column.

`partial = FALSE` refuses a hole rather than returning one: it errors on
a named year that does not reach everywhere, and on ground that no year
has ever covered. Within a year, each tile still takes its latest
flight - one year's flights are spread over several dates between
neighbouring tiles, so years are chosen, not dates.

Orthophotos are the exception to the NoData rule. DOP20 declares no
NoData value, so the uncovered part of an `"rgb"` result comes back
**black** rather than flagged, which no downstream function can tell
from a very dark pixel. Use `partial = FALSE` there, or mask the result
yourself.

## Licence

LGLN open geodata, licensed CC BY 4.0 under section 7 of LGLN's terms of
use. Credit them as **© GeoBasis-DE/LGLN** followed by the year you
downloaded the data, and add **", Daten geändert"** (data modified) for
anything derived from them.

## See also

[`rvt_data_lgln_years()`](https://wiesehahn.github.io/rvtr/reference/rvt_data_lgln_years.md)

## Examples

``` r
# \donttest{
pt <- c(9.9464, 51.6317)                      # Burg Hardenberg
rvt_data_lgln_years(pt, "rgb")
#>   product year first_date  last_date tiles coverage covers
#> 1     rgb 2013 2013-06-18 2013-06-18     1        1   TRUE
#> 2     rgb 2016 2016-03-17 2016-03-17     1        1   TRUE
#> 3     rgb 2019 2019-04-01 2019-04-01     1        1   TRUE
#> 4     rgb 2022 2022-03-13 2022-03-13     1        1   TRUE
#> 5     rgb 2025 2025-03-08 2025-03-08     1        1   TRUE
rvt_data_lgln(pt, "dtm") |> rvt_hillshade() |> rvt_plot()
#> LGLN open data: © GeoBasis-DE/LGLN 2026, CC BY 4.0. Add ", Daten geändert" when you publish anything derived from it.

rvt_data_lgln(pt, "rgb", year = 2013, res = 1) |>
  rvt_plot(range = c(0, 255))

# }
```
