#include <Rcpp.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Geomorphons (Jasiewicz & Stepinski 2013): classify every cell into one of
// ten landform types from the pattern of what it can see around it.
//
// Along each of eight directions the kernel finds the highest and the lowest
// line-of-sight angle within the search radius. Whichever of the two is the
// larger in magnitude decides that direction's ternary element - +1 if the
// terrain rises away, -1 if it falls away, 0 if neither is steep enough to
// clear the flatness threshold. Counting the +1s and -1s over the eight
// directions gives a pair that indexes a lookup table of forms: eight falling
// directions is a peak, eight rising is a pit, and everything between is a
// ridge, spur, slope, hollow, footslope and so on.
//
// Structurally this is the horizon kernel again - same ray walk, same offset
// layout - it just keeps both extremes per direction rather than the maximum
// alone, and never averages them.
//
// Table rows are the count of -1 elements, columns the count of +1, which is
// the convention the GRASS and reference R implementations use.
static const int FORMS10[9][9] = {
  /*        n+ = 0   1   2   3   4   5   6   7   8   */
  /* n- 0 */ {  1,  1,  1,  8,  8,  9,  9,  9, 10 },
  /* n- 1 */ {  1,  1,  8,  8,  8,  9,  9,  9,  0 },
  /* n- 2 */ {  1,  4,  6,  6,  7,  7,  9,  0,  0 },
  /* n- 3 */ {  4,  4,  6,  6,  6,  7,  0,  0,  0 },
  /* n- 4 */ {  4,  4,  5,  6,  6,  0,  0,  0,  0 },
  /* n- 5 */ {  3,  3,  5,  5,  0,  0,  0,  0,  0 },
  /* n- 6 */ {  3,  3,  3,  0,  0,  0,  0,  0,  0 },
  /* n- 7 */ {  3,  3,  0,  0,  0,  0,  0,  0,  0 },
  /* n- 8 */ {  2,  0,  0,  0,  0,  0,  0,  0,  0 }
};

// [[Rcpp::export]]
NumericMatrix geomorphons_kernel(NumericMatrix padded,
                                 int pad,
                                 int nrow_out,
                                 int ncol_out,
                                 IntegerVector dx,
                                 IntegerVector dy,
                                 NumericVector dist,
                                 IntegerVector dir_start,
                                 IntegerVector dir_end,
                                 double flat_threshold,
                                 double flat_distance,
                                 int threads) {

  const int ndir = dir_start.size();
  const int pnrow = padded.nrow();
  const double *P = &padded[0];

  const int *DX = &dx[0];
  const int *DY = &dy[0];
  const double *D = &dist[0];
  const int *DS = &dir_start[0];
  const int *DE = &dir_end[0];

  // Beyond `flat_distance` the fixed angular threshold is replaced by the
  // angle that the same height difference would subtend at that range, so a
  // long search doesn't declare every distant slope significant.
  const double flat_height = flat_distance * std::tan(flat_threshold);

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

      int n_pos = 0, n_neg = 0;

      for (int d = 0; d < ndir; ++d) {
        double zenith = -M_PI / 2.0, nadir = M_PI / 2.0;
        double zdist = 0.0, ndist = 0.0;
        bool any = false;

        const int k1 = DE[d];
        for (int k = DS[d]; k < k1; ++k) {
          const double v = P[(pi + DY[k]) + (R_xlen_t)(pj + DX[k]) * pnrow];
          if (ISNAN(v)) continue;
          any = true;
          const double ang = std::atan2(v - centre, D[k]);
          if (ang > zenith) { zenith = ang; zdist = D[k]; }
          if (ang < nadir)  { nadir  = ang; ndist = D[k]; }
        }
        if (!any) continue;                 // nothing visible: counts as flat

        double zthr = flat_threshold, nthr = flat_threshold;
        if (flat_distance > 0.0) {
          if (zdist > flat_distance) zthr = std::atan2(flat_height, zdist);
          if (ndist > flat_distance) nthr = std::atan2(flat_height, ndist);
        }

        // Significance decides the sign; magnitude only arbitrates when both
        // extremes clear the threshold. Deciding on magnitude alone looks
        // equivalent but is not: on uniformly sloping ground the two angles
        // are exactly equal, the comparison ties, and the direction goes
        // uncounted - which silently turns a cone's summit into a shoulder.
        const bool z_sig = zenith > zthr;
        const bool n_sig = nadir < -nthr;
        if (z_sig && !n_sig) {
          ++n_pos;
        } else if (n_sig && !z_sig) {
          ++n_neg;
        } else if (z_sig && n_sig) {
          const double az = std::fabs(zenith), an = std::fabs(nadir);
          if (az > an) ++n_pos;
          else if (an > az) ++n_neg;
        }
      }

      const int form = FORMS10[n_neg][n_pos];
      OUT[out_idx] = form == 0 ? NA_REAL : (double)form;
    }
  }

  return out;
}
