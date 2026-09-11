#ifndef RVTR_LEVELS_H
#define RVTR_LEVELS_H

#include <Rcpp.h>
#include <cmath>
#include <vector>

// Multi-resolution horizon search, shared by every kernel that needs a
// steepest-angle-per-direction scan: sky-view factor and openness
// (horizon_pyramid.cpp), time in daylight, sky illumination and cast shadow.
//
// Level 0 is the tile itself at full resolution. Each further level is a
// coarser copy of the same DEM covering the distances the previous level
// stopped at, so reaching D map units costs O(log D) instead of O(D). The
// horizon is a max over samples, so splitting a ray across levels is just a
// max over per-level maxima - no reweighting, no double counting. Every angle
// is measured from the cell's own full-resolution elevation; the observer is
// never coarsened.
//
// A coarse level is indexed globally: an output cell at tile offset (i, j) is
// full-resolution pixel (tile_y0 + i, tile_x0 + j), hence coarse cell
// ((tile_y0 + i) / fac, (tile_x0 + j) / fac), less that window's own origin.
// Samples are clamped to the window, which is only ever reached where the
// raster itself ends (see .process_tiled()), so clamping is edge-replication
// at the raster edge and never a tiling artefact.
//
// With no coarse levels at all this reduces to the ordinary single-resolution
// search, so callers need only one code path.

#define RVTR_MAX_LEVELS 16

struct RvtLevels {
  const double *P;            // level 0: the padded tile
  int pnrow, pad;
  int tile_x0, tile_y0;
  int ndir;

  int nlev;                   // coarse levels, beyond level 0
  std::vector<const double*> M;
  std::vector<int> nr, nc, fac, cx0, cy0;

  // offsets per level, level 0 first, so index l + 1 is coarse level l
  std::vector<const int*> DX, DY, DS, DE;
  std::vector<std::vector<double> > IDIST;

  // the R objects above point into; held so they outlive the pointers
  std::vector<Rcpp::IntegerVector> keep_i;
  std::vector<Rcpp::NumericVector> keep_d;
  std::vector<Rcpp::NumericMatrix> keep_m;

  void init(Rcpp::NumericMatrix padded, int pad_, int tx0, int ty0,
            Rcpp::List aux_mats, Rcpp::IntegerVector aux_fac,
            Rcpp::IntegerVector aux_cx0, Rcpp::IntegerVector aux_cy0,
            Rcpp::List off_dx, Rcpp::List off_dy, Rcpp::List off_dist,
            Rcpp::List off_start, Rcpp::List off_end) {
    P = &padded[0];
    pnrow = padded.nrow();
    pad = pad_;
    tile_x0 = tx0;
    tile_y0 = ty0;

    nlev = aux_mats.size();
    if (nlev > RVTR_MAX_LEVELS)
      Rcpp::stop("at most %d coarse levels are supported", RVTR_MAX_LEVELS);

    Rcpp::IntegerVector s0 = off_start[0];
    ndir = s0.size();

    const int nall = nlev + 1;
    DX.resize(nall); DY.resize(nall); DS.resize(nall); DE.resize(nall);
    IDIST.resize(nall);
    keep_i.reserve(nall * 4);
    keep_d.reserve(nall);
    for (int l = 0; l < nall; ++l) {
      Rcpp::IntegerVector dx = off_dx[l], dy = off_dy[l];
      Rcpp::IntegerVector st = off_start[l], en = off_end[l];
      Rcpp::NumericVector di = off_dist[l];
      if (st.size() != ndir)
        Rcpp::stop("every level must supply offsets for the same directions");
      keep_i.push_back(dx); keep_i.push_back(dy);
      keep_i.push_back(st); keep_i.push_back(en);
      keep_d.push_back(di);
      DX[l] = &dx[0]; DY[l] = &dy[0]; DS[l] = &st[0]; DE[l] = &en[0];
      // a reciprocal turns a division per sample per cell into a multiply
      IDIST[l].resize(di.size());
      for (R_xlen_t k = 0; k < di.size(); ++k) IDIST[l][k] = 1.0 / di[k];
    }

    M.resize(nlev); nr.resize(nlev); nc.resize(nlev);
    fac.resize(nlev); cx0.resize(nlev); cy0.resize(nlev);
    keep_m.reserve(nlev);
    for (int l = 0; l < nlev; ++l) {
      Rcpp::NumericMatrix m = aux_mats[l];
      keep_m.push_back(m);
      M[l] = &m[0]; nr[l] = m.nrow(); nc[l] = m.ncol();
      fac[l] = aux_fac[l]; cx0[l] = aux_cx0[l]; cy0[l] = aux_cy0[l];
    }
  }

  // Steepest angle to the skyline in each direction, in radians, measured from
  // `centre`. Writes `ndir` values. Angles may be negative where the terrain
  // falls away; callers that want a half-dome clamp at 0 themselves.
  inline void horizon(int i, int j, double centre, double *out) const {
    const int pi = i + pad, pj = j + pad;

    // where this cell sits in each coarse level - the same for every
    // direction, so resolved once per cell rather than once per ray
    int ci[RVTR_MAX_LEVELS], cj[RVTR_MAX_LEVELS];
    for (int l = 0; l < nlev; ++l) {
      ci[l] = (tile_y0 + i) / fac[l] - cy0[l];
      cj[l] = (tile_x0 + j) / fac[l] - cx0[l];
    }

    for (int d = 0; d < ndir; ++d) {
      double ms = -INFINITY;

      const int k1 = DE[0][d];
      for (int k = DS[0][d]; k < k1; ++k) {
        const double v = P[(pi + DY[0][k]) + (R_xlen_t)(pj + DX[0][k]) * pnrow];
        if (ISNAN(v)) continue;
        const double s = (v - centre) * IDIST[0][k];
        if (s > ms) ms = s;
      }

      for (int l = 0; l < nlev; ++l) {
        const int lv = l + 1, n_r = nr[l], n_c = nc[l];
        const double *MM = M[l];
        const int k2 = DE[lv][d];
        for (int k = DS[lv][d]; k < k2; ++k) {
          int ri = ci[l] + DY[lv][k];
          int rj = cj[l] + DX[lv][k];
          ri = ri < 0 ? 0 : (ri >= n_r ? n_r - 1 : ri);
          rj = rj < 0 ? 0 : (rj >= n_c ? n_c - 1 : rj);
          const double v = MM[ri + (R_xlen_t)rj * n_r];
          if (ISNAN(v)) continue;
          const double s = (v - centre) * IDIST[lv][k];
          if (s > ms) ms = s;
        }
      }

      out[d] = std::atan(ms);
    }
  }
};

#endif
