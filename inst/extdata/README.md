# Reference data

A 1000 x 1000 tile used in examples and as a regression baseline for this
package.

| file | what | committed to git? |
|---|---|---|
| `dtm1.tif` | input DTM, 1 m resolution, ETRS89 / UTM 32N | yes |
| `rvt_sky_view_factor.tif` | rvt-py sky-view factor | no — local dev only |
| `rvt_asvf.tif` | rvt-py anisotropic sky-view factor | no — local dev only |
| `rvt_openness_positive.tif` | rvt-py positive openness | no — local dev only |
| `rvt_openness_negative.tif` | rvt-py negative openness | no — local dev only |
| `rvt_local_dominance.tif` | rvt-py local dominance | no — local dev only |
| `rvt_msrm.tif` | rvt-py multi-scale relief model | no — local dev only |
| `rvt_slrm.tif` | rvt-py simple local relief model | no — local dev only |
| `rvt_slope.tif` | rvt-py slope, degrees | no — local dev only |
| `rvt_hillshade.tif` | rvt-py hillshade, 315°/35° | no — local dev only |
| `rvt_mstp.tif` | rvt-py multi-scale topographic position, 3-band | no — local dev only |
| `rvt_simfilled_dem.tif` | `dtm1.tif` with its NoData holes filled | no — local dev only |
| `rvt_sim_uniform.tif`, `rvt_sim_overcast.tif`, `rvt_sim_uniform_r20.tif` | rvt-py sky illumination | no — local dev only |
| `rvt_shadow.tif` | rvt-py cast shadow, 315°/15° | no — local dev only |
| `rvt_archaeological_vat_combined.tif` | rvt-py archaeological VAT (combined) | no — local dev only |

## Provenance

`dtm1.tif` was rasterised from a Lower Saxony 2017 airborne laser scanning tile
(`3dm_32_557_5699_1_ni_2017.laz`) with `lasR`.

`rvt_sky_view_factor.tif`, `rvt_openness_positive.tif` and
`rvt_archaeological_vat_combined.tif` were produced by the `rvt-py` package
(Relief Visualization Toolbox) at its defaults — 16 directions, 10 px search
radius, no noise removal — via `code/indices/calc_rvt.py` in the sibling
`terrain_visualization` project, then recompressed to `DEFLATE` + `PREDICTOR=3`
when copied here (pixel values unchanged, verified identical not merely close).

The remaining `rvt_*` files were generated directly with `rvt.vis` at its
defaults, via GDAL's Python bindings, uncompressed (not run through
`calc_rvt.py`, which doesn't cover them).

Two exceptions to "at its defaults":

`rvt_asvf.tif` was generated asking rvt-py for `asvf_dir = 45`, not the
default 315. rvt-py attaches its anisotropy weights to the mirrored azimuth,
so its 45 is our 315 — see "Deviations from rvt-py" in `CLAUDE.md`.
Requesting 45 lets the test compare the weighting maths without the deliberate
direction fix getting in the way.

The sky illumination and shadow references were generated from
`rvt_simfilled_dem.tif`, a copy of `dtm1.tif` with its 393 NoData cells
replaced by the mean. rvt-py's `sky_illumination()` pads its DEM pyramid with
`np.min(dem)`, which is `NaN` whenever the DEM has any NoData at all, so it
returns an entirely NaN raster for `dtm1.tif` as-is; filling the holes is the
only way to get a reference out of it. `rvt_sim_uniform_r20.tif` additionally
uses `max_fine_radius = 20`, which equals rvt-py's `max_pyramid_radius` and so
suppresses the DEM pyramid entirely — that is the one configuration where
rvt-py does the same single-resolution search we do, and the two then agree to
a correlation of 0.997.

## Why `dtm1.tif` is here, and the `rvt_*` files aren't

`dtm1.tif` is committed because the package examples (`?rvt_svf` etc.) and the
structural-invariant tests (mosaic identity, tiling identity) in
`tests/testthat/test-reference.R` need something to run against.

The `rvt_*` files are gitignored (`inst/extdata/rvt_*` in `.gitignore`) — this
package doesn't chase bit-for-bit parity with rvt-py (see the main README), so
they're a local development aid for spot-checking, not a shipped fixture. The
"agrees with rvt-py" tests skip themselves when these files are missing. To run
them, regenerate the files with `rvt-py` (see Provenance above) and drop them
back in this folder.

## Alpine tiles (`dgm1_*.tif`) — mountain terrain and a real mosaic

| file | what | committed to git? |
|---|---|---|
| `dgm1_659_5258.tif`, `dgm1_660_5258.tif`, `dgm1_659_5259.tif`, `dgm1_660_5259.tif` | four adjacent 1 km DTM tiles around the Partnachklamm, 1 m, ETRS89 / UTM 32N | no — local dev only |

`dtm1.tif` is gentle lowland: 71 m of relief over a kilometre. That makes it a
poor test of anything where *distant* terrain matters, and this package now has
several such features (`reach` and its multi-resolution search, `rvt_daylight`,
`rvt_shadow`). These four tiles are the opposite case — a gorge in the
Wetterstein with **494 m of relief** over the same footprint.

How much difference that makes, measured as positive openness at a 25 m reach
against a 400 m one, away from the raster edge:

| | mean abs difference | max |
|---|---|---|
| `dtm1.tif` (Lower Saxony) | 0.197° | 2.31° |
| `dgm1_*` mosaic (Partnachklamm) | **6.220°** | **47.56°** |

They also form a genuine four-file 2 × 2 mosaic (2000 × 2000 m, no NoData),
which the mosaic-identity test uses — the `dtm1.tif` version of that test
splits one file into quadrants, so it cannot catch anything that only goes
wrong across real file boundaries.

### Fetching them

Open data from the Bavarian survey authority
(<https://geodaten.bayern.de/opengeodata/>, CC BY 4.0). To restore them:

```r
tiles <- c("659_5259", "660_5259", "659_5258", "660_5258")
for (t in tiles) {
  download.file(paste0("https://download1.bayernwolke.de/a/dgm/dgm1/", t, ".tif"),
                file.path("inst/extdata", paste0("dgm1_", t, ".tif")), mode = "wb")
}
```

Gitignored (`inst/extdata/dgm1_*`) for the same reason as the `rvt_*` files —
about 15 MB of data that is reproducible in one command. The tests that use
them skip themselves when they are missing.

## LGLN sample data — not in this folder at all

The vignette, the relighting examples and the network tests use Lower Saxony
open data fetched with `rvt_data_lgln()`, usually the 1 km tile beside Burg
Hardenberg (EPSG:25832, 565000–566000 E, 5720000–5721000 N):

```r
pt  <- c(9.9464, 51.6317)
dtm <- rvt_data_lgln(pt, "dtm")                        # DGM1, 1 m
dsm <- rvt_data_lgln(pt, "dsm")                        # DOM1, 1 m
rgb <- rvt_data_lgln(pt, "rgb", year = 2025, res = 1)  # DOP20, on the terrain grid
```

Nothing is committed or shipped. By default the result is a virtual raster
over LGLN's remote Cloud-Optimized GeoTIFFs, found through its STAC APIs
(`dgm`, `dom`, `dop` at `<name>.stac.lgln.niedersachsen.de`); `download = TRUE`
caches a copy in `tools::R_user_dir("rvtr", "cache")/lgln/`.

LGLN open geodata, CC BY 4.0 under section 7 of LGLN's terms of use (AGNB).
Credit as `© GeoBasis-DE/LGLN <year of download>`, adding `, Daten geändert`
for anything derived from the data. `tests/testthat/test-data_lgln.R` skips
its network tests when offline.
