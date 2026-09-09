#include <Rcpp.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Local dominance (Kokalj/Hesse, "LiVT"): for each pixel, sum how much an
// observer standing on it (eye height = elevation + observer_height) looks
// down on the terrain sampled at a grid of (distance, angle) offsets, then
// normalise. Unlike the horizon kernel this is a plain accumulation - no
// running max - so every offset contributes independently and the loop
// parallelises the same way (columns independent, no shared state).
//
// [[Rcpp::export]]
NumericMatrix local_dominance_kernel(NumericMatrix padded,
                                     int pad,
                                     int nrow_out,
                                     int ncol_out,
                                     IntegerVector dx,
                                     IntegerVector dy,
                                     NumericVector dist,
                                     NumericVector weight,
                                     double observer_height,
                                     int threads) {

  const int K = dx.size();
  const int pnrow = padded.nrow();
  const double *P = &padded[0];

  const int *DX = &dx[0];
  const int *DY = &dy[0];
  const double *DIST = &dist[0];
  const double *W = &weight[0];

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

      const double eye = centre + observer_height;
      double sum = 0.0;

      for (int k = 0; k < K; ++k) {
        const double v = P[(pi + DY[k]) + (R_xlen_t)(pj + DX[k]) * pnrow];
        if (ISNAN(v)) continue;             // no sample there - contributes nothing
        if (eye > v) sum += (eye - v) / DIST[k] * W[k];
      }

      OUT[out_idx] = sum;
    }
  }

  return out;
}
