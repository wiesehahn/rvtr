#include <Rcpp.h>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Maximum deviation from mean elevation, over a range of window radii.
//
// For one radius, DEV = (z - mean) / sd, both taken over a square window -
// how far the cell sits above or below its surroundings, in units of how
// varied those surroundings are. The kernel evaluates that at every radius
// from `r_min` to `r_max` stepping `r_step`, and keeps whichever result has
// the largest magnitude, so each cell reports the scale at which it stands
// out most.
//
// Mean and sd both come from summed-area tables, so a radius of 2000 costs
// the same as a radius of 3 - which is what makes the broad scale tractable
// at all. Three tables are needed: count and sum (as for MSRM) plus sum of
// squares, since sd^2 = E[z^2] - E[z]^2.
//
// [[Rcpp::export]]
NumericMatrix max_deviation_kernel(NumericMatrix padded,
                                   int pad,
                                   int nrow_out,
                                   int ncol_out,
                                   int r_min,
                                   int r_max,
                                   int r_step,
                                   int threads) {

  const int pnrow = padded.nrow();
  const int pncol = padded.ncol();
  const double *P = &padded[0];

  const R_xlen_t H1 = pnrow + 1, W1 = pncol + 1;
  std::vector<double> sumII((size_t)(H1 * W1), 0.0);
  std::vector<double> sqII((size_t)(H1 * W1), 0.0);
  std::vector<double> cntII((size_t)(H1 * W1), 0.0);

  // Shift elevations towards zero before squaring. Variance doesn't care
  // about a constant offset, but `E[z^2] - E[z]^2` does: on flat ground at
  // 300 m the two terms agree to eight digits before they differ, so the
  // subtraction throws away most of the precision in the answer - and that
  // answer is the standard deviation this metric divides by. Subtracting the
  // mean first keeps the squares small and the cancellation mild. Both terms
  // of the numerator shift equally, so nothing has to be added back.
  double offset = 0.0;
  {
    double s = 0.0;
    R_xlen_t n = 0;
    const R_xlen_t total = (R_xlen_t)pnrow * pncol;
    for (R_xlen_t k = 0; k < total; ++k)
      if (!ISNAN(P[k])) { s += P[k]; ++n; }
    if (n > 0) offset = s / n;
  }

  for (int c = 1; c <= pncol; ++c) {
    for (int r = 1; r <= pnrow; ++r) {
      const double raw = P[(r - 1) + (R_xlen_t)(c - 1) * pnrow];
      const bool ok = !ISNAN(raw);
      const double v = ok ? raw - offset : 0.0;
      const R_xlen_t at = r + (R_xlen_t)c * H1;
      const R_xlen_t up = (r - 1) + (R_xlen_t)c * H1;
      const R_xlen_t lf = r + (R_xlen_t)(c - 1) * H1;
      const R_xlen_t ul = (r - 1) + (R_xlen_t)(c - 1) * H1;
      sumII[at] = v                  + sumII[up] + sumII[lf] - sumII[ul];
      sqII[at]  = v * v              + sqII[up]  + sqII[lf]  - sqII[ul];
      cntII[at] = (ok ? 1.0 : 0.0)   + cntII[up] + cntII[lf] - cntII[ul];
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
      const double centre_raw = P[pi + (R_xlen_t)pj * pnrow];
      const R_xlen_t out_idx = i + (R_xlen_t)j * nrow_out;

      if (ISNAN(centre_raw)) {
        OUT[out_idx] = NA_REAL;
        continue;
      }
      const double centre = centre_raw - offset;   // same shifted space

      double best = 0.0;
      bool have = false;

      for (int r = r_min; r <= r_max; r += r_step) {
        const R_xlen_t r0 = pi - r, r1 = pi + r + 1, c0 = pj - r, c1 = pj + r + 1;
        const double n = cntII[r1 + c1 * H1] - cntII[r0 + c1 * H1]
                       - cntII[r1 + c0 * H1] + cntII[r0 + c0 * H1];
        if (n <= 0) continue;
        const double s = sumII[r1 + c1 * H1] - sumII[r0 + c1 * H1]
                       - sumII[r1 + c0 * H1] + sumII[r0 + c0 * H1];
        const double q = sqII[r1 + c1 * H1] - sqII[r0 + c1 * H1]
                       - sqII[r1 + c0 * H1] + sqII[r0 + c0 * H1];

        const double mean = s / n;
        // fabs guards the tiny negative that rounding can produce on a
        // perfectly flat window; 1e-6 keeps flat ground finite rather than
        // dividing by zero.
        const double sd = std::sqrt(std::fabs(q / n - mean * mean));
        const double dev = (centre - mean) / (sd + 1e-6);

        if (!have || std::fabs(dev) > std::fabs(best)) { best = dev; have = true; }
      }

      OUT[out_idx] = have ? best : NA_REAL;
    }
  }

  return out;
}
