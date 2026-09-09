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
