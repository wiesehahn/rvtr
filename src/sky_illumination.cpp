#include <Rcpp.h>
#include <cmath>
#include <vector>
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
// The search is single-resolution. RVT reaches its 100 px default through a
// DEM pyramid, needed because the equivalent numpy loop is ~10x slower per
// offset than this kernel; here the straight search is affordable and avoids
// approximating the horizon at distance.
//
// [[Rcpp::export]]
NumericMatrix sky_illumination_kernel(NumericMatrix padded,
                                      int pad,
                                      int nrow_out,
                                      int ncol_out,
                                      double xres,
                                      double yres,
                                      IntegerVector dx,
                                      IntegerVector dy,
                                      NumericVector dist,
                                      IntegerVector dir_start,
                                      IntegerVector dir_end,
                                      NumericVector dir_azimuth,
                                      bool overcast,
                                      int threads) {

  const int ndir = dir_start.size();
  const int pnrow = padded.nrow();
  const double *P = &padded[0];

  const int *DX = &dx[0];
  const int *DY = &dy[0];
  const int *DS = &dir_start[0];
  const int *DE = &dir_end[0];

  std::vector<double> inv_dist(dist.size());
  for (R_xlen_t k = 0; k < dist.size(); ++k) inv_dist[k] = 1.0 / dist[k];
  const double *ID = inv_dist.data();

  std::vector<double> az(ndir);
  for (int d = 0; d < ndir; ++d) az[d] = dir_azimuth[d] * M_PI / 180.0;

  // half the angular width of one direction's slice of sky
  const double da = M_PI / (double)ndir;
  const double sin_da = std::sin(da);

  // Flat ground gives exactly this, so dividing by it puts flat ground at 1 -
  // the same convention as sky-view factor and local dominance here. RVT
  // instead divides the overcast model by the brightest pixel in the raster,
  // which cannot survive tiling: the scale would depend on which tile a cell
  // landed in. See the package docs.
  const double flat_norm = overcast
    ? (0.33 * M_PI + 0.67 * (2.0 * M_PI / 3.0))
    : M_PI;

  NumericMatrix out(nrow_out, ncol_out);
  double *OUT = &out[0];

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < ncol_out; ++j) {
    const int pj = j + pad;
    for (int i = 0; i < nrow_out; ++i) {
      const int pi = i + pad;
      const double centre = P[pi + (R_xlen_t)pj * pnrow];
      const R_xlen_t out_idx = i + (R_xlen_t)j * nrow_out;

      if (ISNAN(centre)) {
        OUT[out_idx] = NA_REAL;
        continue;
      }

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

      double sum_a = 0.0, sum_b = 0.0, sum_c = 0.0, sum_d = 0.0;

      for (int d = 0; d < ndir; ++d) {
        // horizon angle, never below the horizontal
        double max_slope = 0.0;
        const int k1 = DE[d];
        for (int k = DS[d]; k < k1; ++k) {
          const double v = P[(pi + DY[k]) + (R_xlen_t)(pj + DX[k]) * pnrow];
          if (ISNAN(v)) continue;
          const double s = (v - centre) * ID[k];
          if (s > max_slope) max_slope = s;
        }
        const double h = std::atan(max_slope);
        const double cos_h = std::cos(h);

        // how much this slice of sky is tipped towards or away from the
        // surface; negative contributions are light the surface cannot see
        const double d_aspect = -2.0 * sin_da * std::cos(az[d] - asp);

        sum_a += cos_h * cos_h;
        const double b = d_aspect * (M_PI / 4.0 - h / 2.0 - std::sin(2.0 * h) / 4.0);
        if (b > 0.0) sum_b += b;

        if (overcast) {
          const double cos3 = cos_h * cos_h * cos_h;
          if (cos3 > 0.0) sum_c += cos3;
          const double dd = d_aspect * (2.0 / 3.0 - cos_h + cos3 / 3.0);
          if (dd > 0.0) sum_d += dd;
        }
      }

      const double cos_slp = std::cos(slp), sin_slp = std::sin(slp);
      const double uniform = da * cos_slp * sum_a
                           + sin_slp * (sum_b < M_PI ? sum_b : M_PI);

      double val;
      if (overcast) {
        // the blend mixes the raw uniform term, not the one already divided
        // by pi - matching the order the reference applies these in
        const double oc = (2.0 * da / 3.0) * cos_slp * sum_c + sin_slp * sum_d;
        val = 0.33 * uniform + 0.67 * oc;
      } else {
        val = uniform;
      }

      OUT[out_idx] = val / flat_norm;
    }
  }

  return out;
}
