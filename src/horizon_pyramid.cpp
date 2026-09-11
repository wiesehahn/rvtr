#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "levels.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Sky-view factor, openness and anisotropic sky-view factor over a
// multi-resolution horizon search. Same outputs as horizon_svf_opns(); the
// scan along each direction is split across levels by RvtLevels (levels.h),
// which carries the reasoning and the indexing rules.
//
// [[Rcpp::export]]
List horizon_pyramid(NumericMatrix padded,
                     int pad,
                     int nrow_out,
                     int ncol_out,
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
                     bool want_svf,
                     bool want_opns,
                     bool want_asvf,
                     NumericVector dir_weight,
                     int threads) {

  RvtLevels L;
  L.init(padded, pad, tile_x0, tile_y0, aux_mats, aux_fac, aux_cx0, aux_cy0,
         off_dx, off_dy, off_dist, off_start, off_end);
  const int ndir = L.ndir;

  NumericMatrix svf(want_svf ? nrow_out : 1, want_svf ? ncol_out : 1);
  NumericMatrix opns(want_opns ? nrow_out : 1, want_opns ? ncol_out : 1);
  NumericMatrix asvf(want_asvf ? nrow_out : 1, want_asvf ? ncol_out : 1);
  double *SVF = &svf[0];
  double *OPNS = &opns[0];
  double *ASVF = &asvf[0];

  // Anisotropy: one weight per direction, largest towards the bright part of
  // the sky. Normalising by their total (rather than by ndir) keeps the
  // result on the same 0-1 scale as plain SVF.
  const double *W = want_asvf ? &dir_weight[0] : NULL;
  double weight_total = 0.0;
  if (want_asvf) for (int d = 0; d < ndir; ++d) weight_total += W[d];

  const double rad2deg = 180.0 / M_PI;

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < ncol_out; ++j) {
    std::vector<double> hor(ndir);      // per thread, reused down the column
    for (int i = 0; i < nrow_out; ++i) {
      const double centre = L.P[(i + pad) + (R_xlen_t)(j + pad) * L.pnrow];
      const R_xlen_t out_idx = i + (R_xlen_t)j * nrow_out;

      if (ISNAN(centre)) {
        if (want_svf) SVF[out_idx] = NA_REAL;
        if (want_opns) OPNS[out_idx] = NA_REAL;
        if (want_asvf) ASVF[out_idx] = NA_REAL;
        continue;
      }

      L.horizon(i, j, centre, hor.data());

      double svf_sum = 0.0, opns_sum = 0.0, asvf_sum = 0.0;
      for (int d = 0; d < ndir; ++d) {
        const double ang = hor[d];
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
