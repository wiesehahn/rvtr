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
// type: 0 profile, 1 plan, 2 tangential, 3 mean, 4 total, 5 gaussian,
//       6 minimal, 7 maximal, 8 unsphericity, 9 casorati, 10 shape index,
//       11 difference, 12 twisting, 13 rotor
//
// Shared by both kernels below - the 3x3 second differences and the fitted
// quadratic differ only in how they arrive at these five derivatives, never in
// what is done with them, so this stays in one place.
static inline double rvt_curv_combine(double zx, double zy,
                                      double zxx, double zyy, double zxy,
                                      int type) {
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

  // The two normal-curvature numerators. Shared by profile, plan and
  // tangential, and by the difference curvature built from the first and
  // third - which is why they are hoisted rather than repeated in the cases.
  const bool flat = (p < 1e-12);
  const double kns_n = zx * zx * zxx + 2.0 * zx * zy * zxy + zy * zy * zyy;
  const double knc_n = zy * zy * zxx - 2.0 * zx * zy * zxy + zx * zx * zyy;
  const double k_prof = flat ? 0.0 : -kns_n / (p * std::pow(q, 1.5));
  const double k_tang = flat ? 0.0 : -knc_n / (p * std::sqrt(q));
  // The twisting numerator, shared by the contour's geodesic torsion and the
  // rotor: the same quantity, normalised two different ways.
  const double tw_n = zx * zy * (zxx - zyy) - zxy * (zx * zx - zy * zy);

  switch (type) {
    case 0: return k_prof;   // profile - along the line of steepest slope
    case 1:                  // plan - the contour through the cell
      return flat ? 0.0 : -knc_n / std::pow(p, 1.5);
    case 2: return k_tang;   // tangential - plan scaled by the slope
    case 3: return H;                            // mean
    // "Total curvature" is a term of art: Evans' un-normalised sum of
    // squared second derivatives, which is also what WhiteboxTools
    // returns. Deliberately *not* k1^2 + k2^2, which would be tidier
    // beside the others here but would quietly mean something different
    // from the same name everywhere else.
    case 4: return zxx * zxx + 2.0 * zxy * zxy + zyy * zyy;   // total
    case 5: return K;                            // gaussian
    case 6: return H - disc;                     // minimal
    case 7: return H + disc;                     // maximal

    // Shary's system again. These four cost nothing: they are arithmetic on
    // H and disc, which every call already computes, or on the two normal
    // curvatures above. Adding them needed no new derivative.
    //
    // unsphericity is (kmax - kmin)/2, i.e. exactly `disc` - how far the
    // surface is from a sphere, zero at an umbilic, never negative. With
    // mean and difference it is one of Shary's three independent components.
    case 8: return disc;
    // Casorati: bending magnitude, sqrt((kmax^2 + kmin^2)/2), which reduces
    // to sqrt(H^2 + disc^2). Zero only on a plane, positive everywhere else,
    // and unlike `total` it is a proper normalised curvature.
    case 9: return std::sqrt(H * H + disc * disc);
    // Koenderink's shape index: where the cell sits on the continuum from
    // pit through saddle to peak, in [-1, 1], independent of how sharply it
    // bends. atan2 rather than atan(H/disc) so the umbilic case (disc = 0)
    // and a plane (both 0) resolve without a guard. Signs follow this
    // package's convex-positive convention, so +1 is a dome and -1 a pit.
    case 10: return (2.0 / M_PI) * std::atan2(H, disc);
    // difference curvature, (kns - knc)/2: which of the two normal
    // curvatures dominates, so it separates flow acceleration from flow
    // convergence where mean curvature merges them.
    case 11: return 0.5 * (k_prof - k_tang);

    // The twisting pair, and the third member of Minar's basic trio
    // (profile, plan, twisting). Both are built from tw_n above; only the
    // normalisation differs, and the normalisation is the whole ballgame -
    // it was verified analytically rather than taken from memory:
    //
    //  - twisting is the contour's geodesic torsion, and equals *exactly*
    //    the rate of change of slope angle (radians per map unit) as you walk
    //    along the contour. Checked against a numerically differentiated
    //    slope angle on z = -x - a.x.y, where it agrees to 8 decimals. That
    //    identity is what pins the (p*q) denominator: q^1.5 or sqrt(q) both
    //    break it.
    //  - rotor is the flow line's curvature seen in plan, hence p^1.5, the
    //    same projected normalisation `plan` uses.
    //
    // A second, independent check: in the principal frame the geodesic
    // torsion is ((k2-k1)/2)*sin(2*theta), so |twisting| can never exceed
    // `unsphericity`. Verified to hold with zero violation over 20000 random
    // derivative draws, and to be attained exactly on a saddle.
    //
    // These two are NOT negated the way profile/plan/mean are. That flip
    // exists to make convex read positive, and twisting has no convex sense
    // to align with - its sign is handedness, not shape. Fixed instead by
    // measurement, on the surface above: positive means "to the right as you
    // look downhill" for both.
    case 12: return flat ? 0.0 : tw_n / (p * q);              // twisting
    default: return flat ? 0.0 : tw_n / std::pow(p, 1.5);     // rotor
  }
}

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

      OUT[out_idx] = rvt_curv_combine(zx, zy, zxx, zyy, zxy, type);
    }
  }

  return out;
}


// Curvature at a chosen scale: the bivariate quadratic
//
//   z = a.dx^2 + b.dy^2 + c.dx.dy + d.dx + e.dy + f
//
// fitted by ordinary least squares over the (2r+1)^2 window, after Evans
// (1980) and Wood (1996). The derivatives come straight off the coefficients
// (zxx = 2a, zyy = 2b, zxy = c, zx = d, zy = e) and then go through exactly
// the same combination as the 3x3 kernel above. dy runs southward, matching
// zy = (down - up) there, so nothing downstream changes sign.
//
// Zevenbergen & Thorne has no n x n form - it interpolates exactly through
// nine points - so a least-squares fit is the only way to widen the window,
// and it is what every multi-scale implementation uses.
//
// Why there is no matrix solve here: the window is always complete, because a
// NoData cell contributes the centre's own elevation (the rule used throughout
// this package). Every odd moment of the design therefore vanishes and the
// normal equations decouple - c, d and e each fall out of a single sum, and
// (a, b, f) reduce to a 2x2 in the sum and difference of a and b. P - Q is
// strictly positive for every r >= 1 by Cauchy-Schwarz, with equality only at
// r = 0, which .cells() already refuses.
//
// At r = 1 this is *not* the kernel above: it averages the second difference
// over all three rows (Evans' method) where that one takes the middle row
// alone. Both are correct; they differ wherever the terrain changes from row
// to row, which is why rvt_curvature() keeps the 3x3 kernel when no radius is
// given and matching software (ArcGIS, WhiteboxTools) uses that form.
//
// Cost is O(r^2) per cell. The six moments are separable - {1, dx, dx^2} along
// rows then {1, dy, dy^2} down columns, the trick log_kernel() already uses
// for the Gaussian - which would make it O(r). Deliberately not done: it needs
// a non-separable fallback wherever the window meets NoData, and at the radii
// this is used at the direct loop is already under the per-output I/O floor.
//
// [[Rcpp::export]]
NumericMatrix curvature_fit_kernel(NumericMatrix padded,
                                   int pad,
                                   int nrow_out,
                                   int ncol_out,
                                   double xres,
                                   double yres,
                                   int r,
                                   int type,
                                   int threads) {

  const int pnrow = padded.nrow();
  const double *P = &padded[0];

  NumericMatrix out(nrow_out, ncol_out);
  double *OUT = &out[0];

  // Design constants: they depend only on the window size, never on the data,
  // so they are built once rather than per cell.
  const double m = 2.0 * r + 1.0;
  double S2 = 0.0, S4 = 0.0;
  for (int k = -r; k <= r; ++k) {
    const double k2 = (double)k * (double)k;
    S2 += k2;
    S4 += k2 * k2;
  }
  const double Pm = m * S4, Qm = S2 * S2, Rm = m * S2, Nm = m * m;
  const double den_sum  = (Pm + Qm) - 2.0 * Rm * Rm / Nm;
  const double den_diff = Pm - Qm;

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < ncol_out; ++j) {
    const int pj = j + pad;
    for (int i = 0; i < nrow_out; ++i) {
      const int pi = i + pad;
      const R_xlen_t out_idx = i + (R_xlen_t)j * nrow_out;
      const double z0 = P[pi + (R_xlen_t)pj * pnrow];

      if (ISNAN(z0)) { OUT[out_idx] = NA_REAL; continue; }

      double Z00 = 0.0, Z10 = 0.0, Z01 = 0.0,
             Z20 = 0.0, Z02 = 0.0, Z11 = 0.0;

      // dx runs east over columns, dy south over rows; the inner loop walks
      // rows, which is the contiguous axis in a column-major matrix
      for (int dx = -r; dx <= r; ++dx) {
        const R_xlen_t col = (R_xlen_t)(pj + dx) * pnrow;
        const double fx = (double)dx;
        for (int dy = -r; dy <= r; ++dy) {
          double v = P[(pi + dy) + col];
          if (ISNAN(v)) v = z0;
          const double fy = (double)dy;
          Z00 += v;
          Z10 += fx * v;
          Z01 += fy * v;
          Z20 += fx * fx * v;
          Z02 += fy * fy * v;
          Z11 += fx * fy * v;
        }
      }

      const double cc = Z11 / Qm;
      const double dd = Z10 / Rm;
      const double ee = Z01 / Rm;
      const double s_ab = (Z20 + Z02 - 2.0 * Rm * Z00 / Nm) / den_sum;
      const double d_ab = (Z20 - Z02) / den_diff;
      const double aa = 0.5 * (s_ab + d_ab);
      const double bb = 0.5 * (s_ab - d_ab);

      OUT[out_idx] = rvt_curv_combine(dd / xres, ee / yres,
                                      2.0 * aa / (xres * xres),
                                      2.0 * bb / (yres * yres),
                                      cc / (xres * yres), type);
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
