#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "levels.h"
#include "skyfactor.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Clear-sky solar irradiation on the terrain, following the ESRA model as
// implemented by GRASS r.sun (Hofierka & Suri 2002; Rigollier et al. 2000).
//
// Two components:
//
//   beam     the direct sun. For each sun position, the cell is lit only if
//            the sun clears its terrain horizon *and* the surface faces the
//            sun at all; the energy that lands is then the direct-normal
//            irradiance times the cosine of the incidence angle. Both of
//            those the daylight kernel already computes - it just throws the
//            magnitudes away and counts.
//
//   diffuse  scattered skylight. Horizontal diffuse irradiance is a function
//            of sun altitude and turbidity alone, so R sums it over the sun
//            positions and passes one number; the only per-cell part is the
//            share of the sky dome this surface can see, which is exactly
//            rvt_sky_factor() - the quantity rvt_sky_illumination() reports.
//            Diffuse does not care whether the *sun* is blocked, so it is not
//            gated by the shadow test.
//
// The attenuation depends on the cell's own elevation through air pressure,
// which is why the beam cannot be fully precomputed in R: at 1200 m the
// atmosphere is some 13% thinner than at sea level, and in mountains that is
// worth several percent of the beam.
//
// [[Rcpp::export]]
List insolation_kernel(NumericMatrix padded,
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
                       NumericVector air_mass0,
                       NumericVector extraterrestrial,
                       NumericVector dir_azimuth,
                       double linke,
                       double diffuse_sum,
                       bool overcast,
                       int threads) {

  RvtLevels L;
  L.init(padded, pad, tile_x0, tile_y0, aux_mats, aux_fac, aux_cx0, aux_cy0,
         off_dx, off_dy, off_dist, off_start, off_end);
  const int ndir = L.ndir;
  const int pnrow = L.pnrow;
  const double *P = L.P;
  const int nsun = sun_altitude.size();

  NumericMatrix beam(nrow_out, ncol_out), diff(nrow_out, ncol_out);
  double *BEAM = &beam[0];
  double *DIFF = &diff[0];

  std::vector<double> az(ndir);
  for (int d = 0; d < ndir; ++d) az[d] = dir_azimuth[d] * M_PI / 180.0;
  const double flat_norm = rvt_sky_flat_norm(overcast);

  // Per sun position, nothing that depends on the cell.
  const double step = 2.0 * M_PI / ndir;
  std::vector<int> bin0(nsun), bin1(nsun);
  std::vector<double> frac(nsun), alt(nsun), sin_alt(nsun), cos_alt(nsun),
                      cos_az(nsun), sin_az(nsun), m0(nsun), g0(nsun);
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
    m0[s] = air_mass0[s];              // relative air mass at sea level
    g0[s] = extraterrestrial[s];       // I0 corrected for Earth-Sun distance
  }

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < ncol_out; ++j) {
    const int pj = j + pad;
    std::vector<double> hor(ndir);

    for (int i = 0; i < nrow_out; ++i) {
      const int pi = i + pad;
      const double centre = P[pi + (R_xlen_t)pj * pnrow];
      const R_xlen_t out_idx = i + (R_xlen_t)j * nrow_out;

      if (ISNAN(centre)) {
        BEAM[out_idx] = NA_REAL;
        DIFF[out_idx] = NA_REAL;
        continue;
      }

      L.horizon(i, j, centre, hor.data());

      double left  = P[pi + (R_xlen_t)(pj - 1) * pnrow];
      double right = P[pi + (R_xlen_t)(pj + 1) * pnrow];
      double up    = P[(pi - 1) + (R_xlen_t)pj * pnrow];
      double down  = P[(pi + 1) + (R_xlen_t)pj * pnrow];
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

      // air pressure relative to sea level at this cell's elevation
      const double pratio = std::exp(-centre / 8434.5);

      double acc = 0.0;
      for (int s = 0; s < nsun; ++s) {
        const double f = frac[s];
        const double h = hor[bin0[s]] * (1.0 - f) + hor[bin1[s]] * f;
        if (alt[s] <= h) continue;                     // behind a hill
        const double cos_rel = cos_asp * cos_az[s] + sin_asp * sin_az[s];
        const double cos_i = sin_alt[s] * cos_slp + cos_alt[s] * sin_slp * cos_rel;
        if (cos_i <= 0.0) continue;                    // facing away

        // ESRA beam: Rayleigh optical thickness at this air mass, then
        // Beer-Lambert attenuation scaled by the Linke turbidity factor
        const double m = pratio * m0[s];
        double dr;
        if (m <= 20.0)
          dr = 1.0 / (6.6296 + m * (1.7513 + m * (-0.1202 + m * (0.0065 - m * 0.00013))));
        else
          dr = 1.0 / (10.4 + 0.718 * m);
        acc += g0[s] * std::exp(-0.8662 * linke * m * dr) * cos_i;
      }

      BEAM[out_idx] = acc;
      DIFF[out_idx] = diffuse_sum *
        rvt_sky_factor(hor.data(), ndir, az.data(), slp, asp, overcast, flat_norm);
    }
  }

  return List::create(_["beam"] = beam, _["diffuse"] = diff);
}
