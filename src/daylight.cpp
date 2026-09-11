#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "levels.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Time in daylight: for how many of the supplied sun positions is each cell
// lit by direct sun?
//
// The efficient formulation is to notice that "lit" is a comparison against a
// horizon that depends only on *azimuth*, never on time:
//
//   lit  <=>  sun_altitude > horizon(sun_azimuth)
//
// where horizon() has two parts, both azimuth-only:
//
//   * the terrain horizon - the steepest angle up to any terrain along that
//     compass bearing, which is what RvtLevels supplies (levels.h), across as
//     many resolutions as the requested reach needs, and
//   * the cell's own surface. A slope facing away from the sun is dark long
//     before the terrain blocks it. Writing the usual incidence test
//     cos(i) = sin(alt)cos(s) + cos(alt)sin(s)cos(az - aspect) > 0 and
//     dividing through by cos(alt) (positive whenever the sun is up) turns it
//     into tan(alt) > -tan(s)cos(az - aspect) - another altitude threshold
//     that varies only with azimuth.
//
// So the expensive part (the terrain search) is done once per cell, and every
// sun position after that is a handful of flops against the stored horizon.
// Adding sun positions - a finer time step, more days - is close to free; it
// is the reach that costs.
//
// Direction d is the compass bearing d * (360/ndir), so the bearing of any sun
// position maps straight onto a bin index with no search. The horizon is
// linearly interpolated between the two bracketing bins rather than snapped to
// the nearest, which keeps low sun angles from banding into visible facets.
//
// [[Rcpp::export]]
NumericMatrix daylight_kernel(NumericMatrix padded,
                              int pad,
                              int nrow_out,
                              int ncol_out,
                              double xres,
                              double yres,
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
                              NumericVector sun_azimuth,
                              NumericVector sun_altitude,
                              int threads) {

  RvtLevels L;
  L.init(padded, pad, tile_x0, tile_y0, aux_mats, aux_fac, aux_cx0, aux_cy0,
         off_dx, off_dy, off_dist, off_start, off_end);
  const int ndir = L.ndir;
  const int nsun = sun_altitude.size();

  if (sun_azimuth.size() != nsun)
    stop("sun_azimuth and sun_altitude must be the same length");

  NumericMatrix out(nrow_out, ncol_out);
  double *OUT = &out[0];

  // Everything about the sun that does not depend on the cell, resolved once.
  const double step = 2.0 * M_PI / ndir;
  std::vector<int> bin0(nsun), bin1(nsun);
  std::vector<double> frac(nsun), alt(nsun), sin_alt(nsun), cos_alt(nsun),
                      cos_az(nsun), sin_az(nsun);
  for (int s = 0; s < nsun; ++s) {
    double a = std::fmod(sun_azimuth[s], 2.0 * M_PI);
    if (a < 0.0) a += 2.0 * M_PI;
    const double b = a / step;
    const int i0 = (int)std::floor(b);
    bin0[s] = i0 % ndir;
    bin1[s] = (i0 + 1) % ndir;
    frac[s] = b - std::floor(b);
    alt[s] = sun_altitude[s];
    sin_alt[s] = std::sin(alt[s]);
    cos_alt[s] = std::cos(alt[s]);
    cos_az[s] = std::cos(a);
    sin_az[s] = std::sin(a);
  }

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < ncol_out; ++j) {
    const int pj = j + pad;
    std::vector<double> hor(ndir);      // per thread, reused down the column

    for (int i = 0; i < nrow_out; ++i) {
      const int pi = i + pad;
      const double centre = L.P[pi + (R_xlen_t)pj * L.pnrow];
      const R_xlen_t out_idx = i + (R_xlen_t)j * nrow_out;

      if (ISNAN(centre)) {
        OUT[out_idx] = NA_REAL;
        continue;
      }

      L.horizon(i, j, centre, hor.data());

      // Local surface, same derivative and edge rule as slope_hillshade().
      double left  = L.P[pi + (R_xlen_t)(pj - 1) * L.pnrow];
      double right = L.P[pi + (R_xlen_t)(pj + 1) * L.pnrow];
      double up    = L.P[(pi - 1) + (R_xlen_t)pj * L.pnrow];
      double down  = L.P[(pi + 1) + (R_xlen_t)pj * L.pnrow];
      if (ISNAN(left))  left  = centre;
      if (ISNAN(right)) right = centre;
      if (ISNAN(up))    up    = centre;
      if (ISNAN(down))  down  = centre;

      double dzdx = ((left - right) / 2.0) / xres;
      double dzdy = ((down - up) / 2.0) / yres;
      if (dzdy == 0.0) dzdy = 10e-9;

      const double slp = std::atan(std::sqrt(dzdx * dzdx + dzdy * dzdy));
      const double asp = std::atan2(dzdx, dzdy);
      const double cos_slp = std::cos(slp), sin_slp = std::sin(slp);
      const double cos_asp = std::cos(asp), sin_asp = std::sin(asp);

      int lit = 0;
      for (int s = 0; s < nsun; ++s) {
        const double f = frac[s];
        const double h = hor[bin0[s]] * (1.0 - f) + hor[bin1[s]] * f;
        if (alt[s] <= h) continue;                    // behind a hill
        // cos(asp - az) expanded, to keep a trig call out of the inner loop
        const double cos_rel = cos_asp * cos_az[s] + sin_asp * sin_az[s];
        if (sin_alt[s] * cos_slp + cos_alt[s] * sin_slp * cos_rel <= 0.0)
          continue;                                   // facing away from the sun
        ++lit;
      }
      OUT[out_idx] = (double)lit;
    }
  }

  return out;
}
