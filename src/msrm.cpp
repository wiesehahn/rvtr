#include <Rcpp.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Multi-scale relief model: differences of box-filtered ("mean") surfaces at
// a sequence of radii, averaged. The box mean at any radius is O(1) per pixel
// via a 2D prefix sum (summed-area table / integral image) built once over
// the tile - this is what makes even a fairly large max radius cheap, rather
// than re-averaging a growing window per radius. A second prefix sum over
// "is this cell valid" makes the mean NoData-aware (NaN cells are excluded,
// not averaged in as zero).
//
// The prefix sums are inherently sequential (each cell depends on its
// neighbours above/left), so that pass is single-threaded - but it is O(N)
// and negligible next to the O(N * num_radii) output pass, which does
// parallelise over columns like the other kernels here.
//
// [[Rcpp::export]]
NumericMatrix msrm_kernel(NumericMatrix padded,
                          int pad,
                          int nrow_out,
                          int ncol_out,
                          IntegerVector radii,
                          int threads) {

  const int pnrow = padded.nrow();
  const int pncol = padded.ncol();
  const double *P = &padded[0];
  const int K = radii.size();
  const int *R = &radii[0];

  const R_xlen_t H1 = pnrow + 1, W1 = pncol + 1;
  std::vector<double> sumII((size_t)(H1 * W1), 0.0);
  std::vector<double> cntII((size_t)(H1 * W1), 0.0);

  for (int c = 1; c <= pncol; ++c) {
    for (int r = 1; r <= pnrow; ++r) {
      const double v = P[(r - 1) + (R_xlen_t)(c - 1) * pnrow];
      const bool ok = !ISNAN(v);
      sumII[r + (R_xlen_t)c * H1] = (ok ? v : 0.0)
        + sumII[(r - 1) + (R_xlen_t)c * H1]
        + sumII[r + (R_xlen_t)(c - 1) * H1]
        - sumII[(r - 1) + (R_xlen_t)(c - 1) * H1];
      cntII[r + (R_xlen_t)c * H1] = (ok ? 1.0 : 0.0)
        + cntII[(r - 1) + (R_xlen_t)c * H1]
        + cntII[r + (R_xlen_t)(c - 1) * H1]
        - cntII[(r - 1) + (R_xlen_t)(c - 1) * H1];
    }
  }

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

      double last_mean = 0.0;
      bool last_valid = false;
      double diff_sum = 0.0;
      int n_diffs = 0;

      for (int k = 0; k < K; ++k) {
        const int r = R[k];
        const R_xlen_t r0 = pi - r, r1 = pi + r + 1, c0 = pj - r, c1 = pj + r + 1;
        const double s = sumII[r1 + c1 * H1] - sumII[r0 + c1 * H1]
                        - sumII[r1 + c0 * H1] + sumII[r0 + c0 * H1];
        const double n = cntII[r1 + c1 * H1] - cntII[r0 + c1 * H1]
                        - cntII[r1 + c0 * H1] + cntII[r0 + c0 * H1];
        const bool valid = n > 0;
        const double mean_k = valid ? s / n : 0.0;

        if (k > 0 && last_valid && valid) {
          diff_sum += (last_mean - mean_k);
          ++n_diffs;
        }
        last_mean = mean_k;
        last_valid = valid;
      }

      OUT[out_idx] = (n_diffs > 0) ? diff_sum / n_diffs : NA_REAL;
    }
  }

  return out;
}
