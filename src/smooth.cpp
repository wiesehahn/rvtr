#include <Rcpp.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Feature-preserving smoothing, after Sun, Rosin, Martin & Langbein (2007),
// "Fast and effective feature-preserving mesh denoising" - the algorithm
// behind WhiteboxTools' FeaturePreservingSmoothing, adapted to a raster.
//
// An ordinary smoothing filter cannot tell noise from a break of slope, so it
// rounds off the terrace edges and gully lips that the metrics here exist to
// find. This works on the surface *normals* instead, in three steps:
//
//   1. a unit normal per cell, from the same central differences the slope
//      kernel uses;
//   2. each normal is replaced by a weighted average of its neighbours',
//      where the weight falls to zero once a neighbour's normal differs by
//      more than `norm_diff`. Across a break of slope the two sides simply
//      do not see each other, so the edge survives while noise - which is
//      small, random deviations within the threshold - averages away;
//   3. elevations are then nudged, repeatedly, towards agreement with those
//      smoothed normals. For a plane through neighbour j with normal n_j,
//      the height at this cell should satisfy n_j . (p_j - p_i) = 0, so each
//      neighbour proposes a correction and the cell takes their mean.
//
// The total elevation change is capped at `max_diff`, which bounds how far
// the result can depart from measured ground.
//
// Every stage reads a bounded neighbourhood, so the whole thing is tileable:
// the caller pads by (iterations + 2) * radius + 1.
//
// [[Rcpp::export]]
NumericMatrix smooth_kernel(NumericMatrix padded,
                            int pad,
                            int nrow_out,
                            int ncol_out,
                            double xres,
                            double yres,
                            int radius,
                            double norm_diff_deg,
                            int iterations,
                            double max_diff,
                            int threads) {

  const int pnrow = padded.nrow();
  const int pncol = padded.ncol();
  const R_xlen_t npad = (R_xlen_t)pnrow * pncol;
  const double *P = &padded[0];

  // Work on a padded copy of the elevations: the update is iterative, and
  // each pass has to read the previous pass everywhere including the border.
  std::vector<double> z(P, P + npad), znew(npad);
  std::vector<double> nx(npad), ny(npad), nz(npad);
  std::vector<double> sx(npad), sy(npad), sz(npad);

  const double cos_thr = std::cos(norm_diff_deg * M_PI / 180.0);

  // ---- 1. unit normals -----------------------------------------------------
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 1; j < pncol - 1; ++j) {
    for (int i = 1; i < pnrow - 1; ++i) {
      const R_xlen_t k = i + (R_xlen_t)j * pnrow;
      const double c = z[k];
      if (ISNAN(c)) { nx[k] = NA_REAL; ny[k] = NA_REAL; nz[k] = NA_REAL; continue; }

      double l = z[k - pnrow], r = z[k + pnrow], u = z[k - 1], d = z[k + 1];
      if (ISNAN(l)) l = c;
      if (ISNAN(r)) r = c;
      if (ISNAN(u)) u = c;
      if (ISNAN(d)) d = c;

      const double dzdx = ((l - r) / 2.0) / xres;
      const double dzdy = ((d - u) / 2.0) / yres;
      const double len = std::sqrt(dzdx * dzdx + dzdy * dzdy + 1.0);
      nx[k] = -dzdx / len; ny[k] = -dzdy / len; nz[k] = 1.0 / len;
    }
  }

  // ---- 2. feature-preserving normal filter ---------------------------------
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = radius + 1; j < pncol - radius - 1; ++j) {
    for (int i = radius + 1; i < pnrow - radius - 1; ++i) {
      const R_xlen_t k = i + (R_xlen_t)j * pnrow;
      if (ISNAN(nz[k])) { sx[k] = NA_REAL; sy[k] = NA_REAL; sz[k] = NA_REAL; continue; }

      double ax = 0.0, ay = 0.0, az = 0.0;
      for (int dj = -radius; dj <= radius; ++dj) {
        for (int di = -radius; di <= radius; ++di) {
          if (di * di + dj * dj > radius * radius) continue;   // circular window
          const R_xlen_t m = k + di + (R_xlen_t)dj * pnrow;
          if (ISNAN(nz[m])) continue;
          const double dot = nx[k] * nx[m] + ny[k] * ny[m] + nz[k] * nz[m];
          if (dot <= cos_thr) continue;        // across a break of slope
          // Sun et al.'s weight: rises smoothly from zero at the threshold,
          // so a neighbour just inside it barely counts
          const double w = (dot - cos_thr) * (dot - cos_thr);
          ax += w * nx[m]; ay += w * ny[m]; az += w * nz[m];
        }
      }
      const double len = std::sqrt(ax * ax + ay * ay + az * az);
      if (len > 0.0) { sx[k] = ax / len; sy[k] = ay / len; sz[k] = az / len; }
      else           { sx[k] = nx[k];    sy[k] = ny[k];    sz[k] = nz[k];    }
    }
  }

  // ---- 3. bring the elevations into line with the smoothed normals ---------
  //
  // Over the *immediate* neighbours only. Sun et al. update a vertex from the
  // faces incident to it, and `radius` belongs to the normal filter alone -
  // letting this step reach the full window instead lets cells on opposite
  // sides of a scarp pull on each other, which smears the very edge the
  // normal filter went to such trouble to keep (measured: a one-cell step
  // spread over eight).
  for (int it = 0; it < iterations; ++it) {
    std::copy(z.begin(), z.end(), znew.begin());
#ifdef _OPENMP
    #pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (int j = radius + 1; j < pncol - radius - 1; ++j) {
      for (int i = radius + 1; i < pnrow - radius - 1; ++i) {
        const R_xlen_t k = i + (R_xlen_t)j * pnrow;
        if (ISNAN(z[k]) || ISNAN(sz[k])) continue;

        double acc = 0.0;
        int cnt = 0;
        for (int dj = -1; dj <= 1; ++dj) {
          for (int di = -1; di <= 1; ++di) {
            const R_xlen_t m = k + di + (R_xlen_t)dj * pnrow;
            if (ISNAN(z[m]) || ISNAN(sz[m])) continue;
            // the plane through neighbour m wants n . (p_m - p_i) = 0; solving
            // for this cell's height gives the correction it proposes
            const double dx = di * xres, dy = dj * yres;
            acc += sz[m] * (sx[m] * dx + sy[m] * dy + sz[m] * (z[m] - z[k]));
            ++cnt;
          }
        }
        if (cnt > 0) znew[k] = z[k] + acc / cnt;
      }
    }
    std::swap(z, znew);
  }

  // ---- crop back, capping how far any cell may have moved ------------------
  NumericMatrix out(nrow_out, ncol_out);
  double *OUT = &out[0];
  for (int j = 0; j < ncol_out; ++j) {
    for (int i = 0; i < nrow_out; ++i) {
      const R_xlen_t k = (i + pad) + (R_xlen_t)(j + pad) * pnrow;
      const double orig = P[k];
      if (ISNAN(orig)) { OUT[i + (R_xlen_t)j * nrow_out] = NA_REAL; continue; }
      double v = z[k];
      if (ISNAN(v)) v = orig;
      const double d = v - orig;
      if (d > max_diff) v = orig + max_diff;
      else if (d < -max_diff) v = orig - max_diff;
      OUT[i + (R_xlen_t)j * nrow_out] = v;
    }
  }
  return out;
}
