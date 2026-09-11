#include <Rcpp.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Curvature, from the Zevenbergen & Thorne (1987) quadratic fitted to the
// 3x3 neighbourhood. The partial derivatives fall straight out of the fit:
//
//   zxx = (left + right)/2 - centre        (all over the cell size squared)
//   zyy = (up + down)/2 - centre
//   zxy = (ne - nw - se + sw) / 4
//   zx  = (right - left) / 2
//   zy  = (down - up) / 2
//
// and the named curvatures are the standard combinations of those. Sign
// convention: positive is convex (a ridge or bank), negative concave (a
// hollow or ditch), which is the usual convention in geomorphometry and the
// one that makes profile curvature read the same way as the relief models
// here.
//
// type: 0 profile, 1 plan, 2 tangential, 3 mean, 4 total, 5 gaussian
//
// [[Rcpp::export]]
NumericMatrix curvature_kernel(NumericMatrix padded,
                               int pad,
                               int nrow_out,
                               int ncol_out,
                               double xres,
                               double yres,
                               int type,
                               int threads) {

  const int pnrow = padded.nrow();
  const double *P = &padded[0];

  NumericMatrix out(nrow_out, ncol_out);
  double *OUT = &out[0];

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < ncol_out; ++j) {
    const int pj = j + pad;
    for (int i = 0; i < nrow_out; ++i) {
      const int pi = i + pad;
      const R_xlen_t out_idx = i + (R_xlen_t)j * nrow_out;
      const double z = P[pi + (R_xlen_t)pj * pnrow];

      if (ISNAN(z)) { OUT[out_idx] = NA_REAL; continue; }

      // NoData neighbours fall back to the centre, as elsewhere in this
      // package, so a hole does not eat a ring of cells around itself
      double zl = P[pi + (R_xlen_t)(pj - 1) * pnrow];
      double zr = P[pi + (R_xlen_t)(pj + 1) * pnrow];
      double zu = P[(pi - 1) + (R_xlen_t)pj * pnrow];
      double zd = P[(pi + 1) + (R_xlen_t)pj * pnrow];
      double znw = P[(pi - 1) + (R_xlen_t)(pj - 1) * pnrow];
      double zne = P[(pi - 1) + (R_xlen_t)(pj + 1) * pnrow];
      double zsw = P[(pi + 1) + (R_xlen_t)(pj - 1) * pnrow];
      double zse = P[(pi + 1) + (R_xlen_t)(pj + 1) * pnrow];
      if (ISNAN(zl)) zl = z;   if (ISNAN(zr)) zr = z;
      if (ISNAN(zu)) zu = z;   if (ISNAN(zd)) zd = z;
      if (ISNAN(znw)) znw = z; if (ISNAN(zne)) zne = z;
      if (ISNAN(zsw)) zsw = z; if (ISNAN(zse)) zse = z;

      // Second differences, not the Zevenbergen & Thorne polynomial
      // coefficients: their D and E are half the corresponding partial
      // derivative, which is why their published curvature formulae carry a
      // leading 2. Using the derivatives directly keeps the expressions below
      // the textbook ones.
      const double zx  = (zr - zl) / (2.0 * xres);
      const double zy  = (zd - zu) / (2.0 * yres);
      const double zxx = (zl + zr - 2.0 * z) / (xres * xres);
      const double zyy = (zu + zd - 2.0 * z) / (yres * yres);
      // y runs southward here, matching zy above, so the "south" corners are
      // the +y ones
      const double zxy = (zse - zsw - zne + znw) / (4.0 * xres * yres);

      const double p = zx * zx + zy * zy;      // squared gradient magnitude
      const double q = p + 1.0;

      // Mean and Gaussian curvature in their proper form - with the (1+p)
      // denominators, not just the numerators. Dropping those is only correct
      // where the ground is level, and silently wrong everywhere else: on a
      // 17-degree slope the bare Laplacian overstates mean curvature by more
      // than a factor of two.
      //
      // Negated so that positive means convex, matching profile/plan below.
      // The textbook H is positive for a bowl; here a bowl is concave and
      // reads negative, which is the convention that makes these agree with
      // the relief models elsewhere in the package.
      const double H = -((1.0 + zy * zy) * zxx - 2.0 * zx * zy * zxy
                         + (1.0 + zx * zx) * zyy) / (2.0 * std::pow(q, 1.5));
      const double K = (zxx * zyy - zxy * zxy) / (q * q);
      // principal curvatures; the radicand is zero at an umbilic and can go
      // slightly negative through rounding
      const double disc = std::sqrt(std::max(0.0, H * H - K));

      double v;
      switch (type) {
        case 0:  // profile - curvature along the line of steepest slope
          v = (p < 1e-12) ? 0.0
            : -(zx * zx * zxx + 2.0 * zx * zy * zxy + zy * zy * zyy) /
               (p * std::pow(q, 1.5));
          break;
        case 1:  // plan - curvature of the contour through the cell
          v = (p < 1e-12) ? 0.0
            : -(zy * zy * zxx - 2.0 * zx * zy * zxy + zx * zx * zyy) /
               std::pow(p, 1.5);
          break;
        case 2:  // tangential - plan curvature scaled by the slope
          v = (p < 1e-12) ? 0.0
            : -(zy * zy * zxx - 2.0 * zx * zy * zxy + zx * zx * zyy) /
               (p * std::sqrt(q));
          break;
        case 3: v = H; break;                       // mean
        // "Total curvature" is a term of art: Evans' un-normalised sum of
        // squared second derivatives, which is also what WhiteboxTools
        // returns. Deliberately *not* k1^2 + k2^2, which would be tidier
        // beside the others here but would quietly mean something different
        // from the same name everywhere else.
        case 4: v = zxx * zxx + 2.0 * zxy * zxy + zyy * zyy; break;
        case 5: v = K; break;                       // gaussian
        case 6: v = H - disc; break;                // minimal
        default: v = H + disc;                      // maximal
      }
      OUT[out_idx] = v;
    }
  }

  return out;
}


// Laplacian of Gaussian: smooth, then take the second derivative. Written as
// two passes rather than one convolution because the Gaussian is separable -
// a horizontal pass then a vertical one - which makes the cost linear in the
// kernel radius instead of quadratic.
//
// [[Rcpp::export]]
NumericMatrix log_kernel(NumericMatrix padded,
                         int pad,
                         int nrow_out,
                         int ncol_out,
                         NumericVector gauss,
                         int threads) {

  const int pnrow = padded.nrow();
  const int pncol = padded.ncol();
  const double *P = &padded[0];
  const int g = gauss.size();
  const int gr = (g - 1) / 2;                  // gaussian kernel radius
  const double *G = &gauss[0];

  // horizontal pass over the whole padded tile, then vertical; the blurred
  // result keeps a one-cell margin for the Laplacian that follows
  NumericMatrix tmp(pnrow, pncol), blur(pnrow, pncol);
  double *T = &tmp[0], *B = &blur[0];

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int c = 0; c < pncol; ++c) {
    for (int r = 0; r < pnrow; ++r) {
      double s = 0.0, w = 0.0;
      for (int k = -gr; k <= gr; ++k) {
        const int cc = c + k;
        if (cc < 0 || cc >= pncol) continue;
        const double v = P[r + (R_xlen_t)cc * pnrow];
        if (ISNAN(v)) continue;
        s += v * G[k + gr];
        w += G[k + gr];
      }
      T[r + (R_xlen_t)c * pnrow] = w > 0.0 ? s / w : NA_REAL;
    }
  }

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int c = 0; c < pncol; ++c) {
    for (int r = 0; r < pnrow; ++r) {
      double s = 0.0, w = 0.0;
      for (int k = -gr; k <= gr; ++k) {
        const int rr = r + k;
        if (rr < 0 || rr >= pnrow) continue;
        const double v = T[rr + (R_xlen_t)c * pnrow];
        if (ISNAN(v)) continue;
        s += v * G[k + gr];
        w += G[k + gr];
      }
      B[r + (R_xlen_t)c * pnrow] = w > 0.0 ? s / w : NA_REAL;
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
      const R_xlen_t out_idx = i + (R_xlen_t)j * nrow_out;
      const double c0 = B[pi + (R_xlen_t)pj * pnrow];
      if (ISNAN(c0) || ISNAN(P[pi + (R_xlen_t)pj * pnrow])) {
        OUT[out_idx] = NA_REAL;
        continue;
      }
      const double l = B[pi + (R_xlen_t)(pj - 1) * pnrow];
      const double r = B[pi + (R_xlen_t)(pj + 1) * pnrow];
      const double u = B[(pi - 1) + (R_xlen_t)pj * pnrow];
      const double d = B[(pi + 1) + (R_xlen_t)pj * pnrow];
      OUT[out_idx] = (ISNAN(l) || ISNAN(r) || ISNAN(u) || ISNAN(d))
        ? NA_REAL : (l + r + u + d - 4.0 * c0);
    }
  }

  return out;
}
