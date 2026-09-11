#ifndef RVTR_SKYFACTOR_H
#define RVTR_SKYFACTOR_H

#include <Rcpp.h>
#include <cmath>

// The share of diffuse sky radiation a cell receives, given its horizon in
// each direction and the way the surface is tilted. Normalised so flat, open
// ground is exactly 1 - so it doubles as the transposition factor that turns
// horizontal diffuse irradiance into irradiance on this particular cell.
//
// Shared by sky_illumination_kernel() and the insolation kernel, which need
// the same integral for different purposes: the first reports it directly, the
// second multiplies horizontal diffuse irradiance by it.
//
// `overcast` switches between a uniformly bright sky (the isotropic
// assumption, right for a clear-sky diffuse transposition) and the standard
// cloudy-day distribution that is brighter overhead.

// Value flat ground produces, which is what the result is divided by.
inline double rvt_sky_flat_norm(bool overcast) {
  return overcast ? (0.33 * M_PI + 0.67 * (2.0 * M_PI / 3.0)) : M_PI;
}

inline double rvt_sky_factor(const double *hor, int ndir, const double *az,
                             double slp, double asp, bool overcast,
                             double flat_norm) {
  const double da = M_PI / (double)ndir;   // half-width of one slice of sky
  const double sin_da = std::sin(da);

  double sum_a = 0.0, sum_b = 0.0, sum_c = 0.0, sum_d = 0.0;
  for (int d = 0; d < ndir; ++d) {
    // horizon angle, never below the horizontal: a cell on a summit sees a
    // full half-dome and no more
    const double h = hor[d] > 0.0 ? hor[d] : 0.0;
    const double cos_h = std::cos(h);

    // how much this slice of sky is tipped towards or away from the surface;
    // negative contributions are light the surface cannot see
    const double d_aspect = -2.0 * sin_da * std::cos(az[d] - asp);

    sum_a += cos_h * cos_h;
    const double b = d_aspect * (M_PI / 4.0 - h / 2.0 - std::sin(2.0 * h) / 4.0);
    if (b > 0.0) sum_b += b;

    if (overcast) {
      const double cos3 = cos_h * cos_h * cos_h;
      if (cos3 > 0.0) sum_c += cos3;
      const double dd = d_aspect * (2.0 / 3.0 - cos_h + cos3 / 3.0);
      if (dd > 0.0) sum_d += dd;
    }
  }

  const double cos_slp = std::cos(slp), sin_slp = std::sin(slp);
  const double uniform = da * cos_slp * sum_a
                       + sin_slp * (sum_b < M_PI ? sum_b : M_PI);

  double val;
  if (overcast) {
    // the blend mixes the raw uniform term, not the one already divided by
    // pi - matching the order the reference applies these in
    const double oc = (2.0 * da / 3.0) * cos_slp * sum_c + sin_slp * sum_d;
    val = 0.33 * uniform + 0.67 * oc;
  } else {
    val = uniform;
  }
  return val / flat_norm;
}

#endif
