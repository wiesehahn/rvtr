#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "levels.h"
#include "skyfactor.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Diffuse sky illumination: how much light a cell receives from the whole sky
// dome, given both its own tilt and the terrain blocking the view around it.
//
// Structurally this is the horizon kernel again - a max slope per direction -
// but each direction's contribution is weighted by where that part of the sky
// sits relative to the surface normal, so slope and aspect enter as well. Two
// sky models are supported: "uniform" (equally bright everywhere) and
// "overcast" (brighter overhead, the standard cloudy-day distribution).
//
// The horizon search itself comes from RvtLevels (levels.h), so a long reach
// is served from progressively coarser copies of the DEM rather than scanned
// at full resolution throughout. RVT reaches its 100 m default through a DEM
// pyramid too, though it needs one because the equivalent numpy loop is ~10x
// slower per offset than this kernel.
//
// [[Rcpp::export]]
NumericMatrix sky_illumination_kernel(NumericMatrix padded,
                                      int pad,
                                      int nrow_out,
                                      int ncol_out,
                                      double xres,
                                      double yres,
                                      int tile_x0,
                                      int tile_y0,
                                      List aux_mats,
                                      IntegerVector aux_fac,
                                      IntegerVector aux_cx0,
                                      IntegerVector aux_cy0,
                                      List off_dx,
                                      List off_dy,
                                      List off_dist,
                                      List off_start,
                                      List off_end,
                                      NumericVector dir_azimuth,
                                      bool overcast,
                                      int threads) {

  RvtLevels L;
  L.init(padded, pad, tile_x0, tile_y0, aux_mats, aux_fac, aux_cx0, aux_cy0,
         off_dx, off_dy, off_dist, off_start, off_end);
  const int ndir = L.ndir;
  const int pnrow = L.pnrow;
  const double *P = L.P;

  std::vector<double> az(ndir);
  for (int d = 0; d < ndir; ++d) az[d] = dir_azimuth[d] * M_PI / 180.0;

  // Flat ground gives exactly this, so dividing by it puts flat ground at 1 -
  // the same convention as sky-view factor and local dominance here. RVT
  // instead divides the overcast model by the brightest pixel in the raster,
  // which cannot survive tiling: the scale would depend on which tile a cell
  // landed in. See the package docs.
  const double flat_norm = rvt_sky_flat_norm(overcast);

  NumericMatrix out(nrow_out, ncol_out);
  double *OUT = &out[0];

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < ncol_out; ++j) {
    const int pj = j + pad;
    std::vector<double> hor(ndir);      // per thread, reused down the column
    for (int i = 0; i < nrow_out; ++i) {
      const int pi = i + pad;
      const double centre = P[pi + (R_xlen_t)pj * pnrow];
      const R_xlen_t out_idx = i + (R_xlen_t)j * nrow_out;

      if (ISNAN(centre)) {
        OUT[out_idx] = NA_REAL;
        continue;
      }

      L.horizon(i, j, centre, hor.data());

      // slope and aspect, same derivative and NoData fallback as the
      // slope/hillshade kernel
      double left  = P[pi + (R_xlen_t)(pj - 1) * pnrow];
      double right = P[pi + (R_xlen_t)(pj + 1) * pnrow];
      double up    = P[(pi - 1) + (R_xlen_t)pj * pnrow];
      double down  = P[(pi + 1) + (R_xlen_t)pj * pnrow];
      if (ISNAN(left))  left  = centre;
      if (ISNAN(right)) right = centre;
      if (ISNAN(up))    up    = centre;
      if (ISNAN(down))  down  = centre;

      double dzdx = ((left - right) / 2.0) / xres;
      double dzdy = ((down - up) / 2.0) / yres;
      if (dzdy == 0.0) dzdy = 10e-9;
      const double slp = std::atan(std::sqrt(dzdx * dzdx + dzdy * dzdy));
      const double asp = std::atan2(dzdx, dzdy);

      OUT[out_idx] = rvt_sky_factor(hor.data(), ndir, az.data(), slp, asp,
                                     overcast, flat_norm);
    }
  }

  return out;
}
