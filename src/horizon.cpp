#include <Rcpp.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Fused horizon search producing sky-view factor and positive openness in a
// single pass over the output, with no intermediate rasters.
//
// The obvious vectorised implementation allocates a full-raster temporary per
// shift step (shift, subtract/divide, running max) - at 4000x4000 that is
// ~150 steps x 3 x 122 MB of memory traffic. Here every direction and every
// offset is resolved per pixel in registers, so the only memory written is
// the two outputs. Reads stay within a +/- radius_max band of rows, which
// stays in cache.
//
// Columns are independent, so the outer loop parallelises with no locking.
//
// [[Rcpp::export]]
List horizon_svf_opns(NumericMatrix padded,
                      int pad,
                      int nrow_out,
                      int ncol_out,
                      IntegerVector dx,
                      IntegerVector dy,
                      NumericVector dist,
                      IntegerVector dir_start,
                      IntegerVector dir_end,
                      bool want_svf,
                      bool want_opns,
                      bool want_asvf,
                      NumericVector dir_weight,
                      int threads) {

  const int ndir = dir_start.size();
  const int pnrow = padded.nrow();
  const double *P = &padded[0];

  NumericMatrix svf(want_svf ? nrow_out : 1, want_svf ? ncol_out : 1);
  NumericMatrix opns(want_opns ? nrow_out : 1, want_opns ? ncol_out : 1);
  NumericMatrix asvf(want_asvf ? nrow_out : 1, want_asvf ? ncol_out : 1);
  double *SVF = &svf[0];
  double *OPNS = &opns[0];
  double *ASVF = &asvf[0];

  const int *DX = &dx[0];
  const int *DY = &dy[0];
  const int *DS = &dir_start[0];
  const int *DE = &dir_end[0];

  // Anisotropy: one weight per direction, largest towards the bright part of
  // the sky. Normalising by their total (rather than by ndir) keeps the
  // result on the same 0-1 scale as plain SVF.
  const double *W = want_asvf ? &dir_weight[0] : NULL;
  double weight_total = 0.0;
  if (want_asvf) for (int d = 0; d < ndir; ++d) weight_total += W[d];

  // reciprocals turn a division per offset per pixel into a multiply
  std::vector<double> inv_dist(dist.size());
  for (R_xlen_t k = 0; k < dist.size(); ++k) inv_dist[k] = 1.0 / dist[k];
  const double *ID = inv_dist.data();

  const double rad2deg = 180.0 / M_PI;

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
        if (want_svf) SVF[out_idx] = NA_REAL;
        if (want_opns) OPNS[out_idx] = NA_REAL;
        if (want_asvf) ASVF[out_idx] = NA_REAL;
        continue;
      }

      double svf_sum = 0.0;
      double opns_sum = 0.0;
      double asvf_sum = 0.0;

      for (int d = 0; d < ndir; ++d) {
        // No valid neighbour in this direction (only reachable beside
        // NoData): treat it as fully open, i.e. the horizon angle is exactly
        // -90 deg. SVF is unaffected either way, since the angle is clamped
        // at 0 before the sine.
        double max_slope = -INFINITY;
        const int k1 = DE[d];
        for (int k = DS[d]; k < k1; ++k) {
          const double v = P[(pi + DY[k]) + (R_xlen_t)(pj + DX[k]) * pnrow];
          if (ISNAN(v)) continue;            // matches numpy fmax semantics
          const double s = (v - centre) * ID[k];
          if (s > max_slope) max_slope = s;
        }
        const double ang = std::atan(max_slope);
        if (want_svf || want_asvf) {
          const double vis = 1.0 - std::sin(ang > 0.0 ? ang : 0.0);
          if (want_svf) svf_sum += vis;
          if (want_asvf) asvf_sum += vis * W[d];
        }
        if (want_opns) opns_sum += ang;
      }

      if (want_svf) SVF[out_idx] = svf_sum / ndir;
      if (want_opns) OPNS[out_idx] = 90.0 - (opns_sum / ndir) * rad2deg;
      if (want_asvf) ASVF[out_idx] = asvf_sum / weight_total;
    }
  }

  return List::create(_["svf"] = want_svf ? svf : NumericMatrix(0, 0),
                      _["opns"] = want_opns ? opns : NumericMatrix(0, 0),
                      _["asvf"] = want_asvf ? asvf : NumericMatrix(0, 0));
}


// Slope (radians) plus one hillshade per requested sun position, sharing the
// same derivative computation. `sun_azimuths` and `sun_elevations` are paired
// element by element, so this covers one sun (rvt_hillshade), several
// elevations from one direction (rvt_vat) and several directions at one
// elevation (rvt_multi_hillshade) without recomputing the derivatives each
// time. Pass empty vectors for slope only.
//
// [[Rcpp::export]]
List slope_hillshade(NumericMatrix padded,
                     int pad,
                     int nrow_out,
                     int ncol_out,
                     double xres,
                     double yres,
                     NumericVector sun_azimuths,
                     NumericVector sun_elevations,
                     int threads) {

  const int pnrow = padded.nrow();
  const double *P = &padded[0];
  const int n_sun = sun_elevations.size();
  if (sun_azimuths.size() != n_sun)
    stop("sun_azimuths and sun_elevations must be the same length");

  NumericMatrix slope(nrow_out, ncol_out);
  double *SLP = &slope[0];

  List hs(n_sun);
  std::vector<double*> HS(n_sun);
  for (int e = 0; e < n_sun; ++e) {
    NumericMatrix m(nrow_out, ncol_out);
    hs[e] = m;
    HS[e] = &(as<NumericMatrix>(hs[e]))[0];
  }

  std::vector<double> cos_z(n_sun), sin_z(n_sun), az_rad(n_sun);
  for (int e = 0; e < n_sun; ++e) {
    const double zenith = M_PI / 2.0 - sun_elevations[e] * M_PI / 180.0;
    cos_z[e] = std::cos(zenith);
    sin_z[e] = std::sin(zenith);
    az_rad[e] = sun_azimuths[e] * M_PI / 180.0;
  }

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < ncol_out; ++j) {
    const int pj = j + pad;
    for (int i = 0; i < nrow_out; ++i) {
      const int pi = i + pad;
      const R_xlen_t out_idx = i + (R_xlen_t)j * nrow_out;

      const double centre = P[pi + (R_xlen_t)pj * pnrow];
      if (ISNAN(centre)) {
        SLP[out_idx] = NA_REAL;
        for (int e = 0; e < n_sun; ++e) HS[e][out_idx] = NA_REAL;
        continue;
      }

      // A NoData neighbour falls back to the centre cell's own elevation,
      // which is the same edge-replicate rule already applied at the raster
      // boundary - just extended to holes in the interior. Without it a valid
      // cell touching a hole would go NoData itself, eroding a one-pixel ring
      // around every hole, and none of the other metrics here do that.
      double left  = P[pi + (R_xlen_t)(pj - 1) * pnrow];
      double right = P[pi + (R_xlen_t)(pj + 1) * pnrow];
      double up    = P[(pi - 1) + (R_xlen_t)pj * pnrow];
      double down  = P[(pi + 1) + (R_xlen_t)pj * pnrow];
      if (ISNAN(left))  left  = centre;
      if (ISNAN(right)) right = centre;
      if (ISNAN(up))    up    = centre;
      if (ISNAN(down))  down  = centre;

      // rvt.vis.slope_aspect(): dzdx from (left - right), dzdy from (below - above)
      double dzdx = ((left - right) / 2.0) / xres;
      double dzdy = ((down - up) / 2.0) / yres;
      if (dzdy == 0.0) dzdy = 10e-9;         // numeric stability, as in RVT

      const double tan_slope = std::sqrt(dzdx * dzdx + dzdy * dzdy);
      const double slp = std::atan(tan_slope);
      const double asp = std::atan2(dzdx, dzdy);

      SLP[out_idx] = slp;
      const double cos_slp = std::cos(slp), sin_slp = std::sin(slp);
      for (int e = 0; e < n_sun; ++e) {
        double v = cos_z[e] * cos_slp
                 + sin_z[e] * sin_slp * std::cos(asp - az_rad[e]);
        HS[e][out_idx] = v < 0.0 ? 0.0 : v;
      }
    }
  }

  return List::create(_["slope"] = slope, _["hillshade"] = hs);
}
