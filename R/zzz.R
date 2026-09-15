## Package load hooks.

.rvtr_env <- new.env(parent = emptyenv())

## Many European survey GeoTIFFs - the Lower Saxony data behind dtm1.tif and
## rvt_data_lgln(), for a start - describe ETRS89 through GeoTIFF keys that
## differ trivially from the EPSG registry's definition. GDAL then warns on
## every open, and because rasters are reopened for every tile, a single call
## floods the console. GTIFF_SRS_SOURCE=EPSG makes GDAL take the registry
## definition, which is the fix GDAL's own message recommends. Measured: the
## warnings go from 6 to 0 on a small script, and the CRS written to outputs,
## and the values computed, come out identical either way.
##
## Only set when the user has not chosen a value themselves (as an option or an
## environment variable), and cleared again on unload so the setting does not
## outlive the package.
.onLoad <- function(libname, pkgname) {
  if (!nzchar(gdalraster::get_config_option("GTIFF_SRS_SOURCE"))) {
    gdalraster::set_config_option("GTIFF_SRS_SOURCE", "EPSG")
    .rvtr_env$set_srs_source <- TRUE
  }
}

.onUnload <- function(libpath) {
  if (isTRUE(.rvtr_env$set_srs_source))
    gdalraster::set_config_option("GTIFF_SRS_SOURCE", "")
}
