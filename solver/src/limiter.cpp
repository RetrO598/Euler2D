#include <pre/parameter.h>
#include <solver/limiter.h>
#include <solver/variableDef.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace solver {
BaseLimiter::BaseLimiter(const preprocess::parameter &param,
                         const preprocess::Geometry &geom,
                         std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
                         std::vector<PRIM_VAR> &umin,
                         std::vector<PRIM_VAR> &umax,
                         std::vector<PRIM_VAR> &lim,
                         std::vector<PRIM_VAR> &gradx,
                         std::vector<PRIM_VAR> &grady)
    : param(param),
      geom(geom),
      cv(cv),
      dv(dv),
      umin(umin),
      umax(umax),
      lim(lim),
      gradx(gradx),
      grady(grady) {
  limRef.dens = 0.0;
  limRef.press = 0.0;
  limRef.velx = 0.0;
  limRef.vely = 0.0;
}

void BaseLimiter::limiterRefVals() {
  volRef = -1.0e+32;
  for (int i = 0; i < geom.phyNodes; ++i) {
    volRef = std::max(volRef, geom.vol[i]);
  }

  double gam1 = param.gamma - 1.0;
  double rgas = gam1 * param.Cp / param.gamma;

  if (param.flowtype_ == preprocess::flowType::External) {
    limRef.dens = param.RhoInfinity;
    limRef.velx = std::sqrt(param.uInfinity * param.uInfinity +
                            param.vInfinity * param.vInfinity);
    limRef.vely = limRef.velx;
    limRef.press = param.PsInfinity;
  } else {
    double temp = param.TtInlet *
                  std::pow(param.PsOutlet / param.PtInlet, gam1 / param.gamma);
    double rho = param.PsOutlet / (rgas * temp);
    double cs = std::sqrt(param.gamma * param.PsOutlet / rho);
    double mach = std::sqrt(2.0 * ((param.TtInlet / temp) - 1.0) / gam1);
    limRef.dens = rho;
    limRef.velx = mach * cs;
    limRef.vely = limRef.velx;
    limRef.press = param.PsOutlet;
  }
}

void BaseLimiter::limiterInit() {
  for (int i = 0; i < geom.phyNodes; ++i) {
    umin[i].dens = cv[i].dens;
    umin[i].velx = cv[i].xmom / cv[i].dens;
    umin[i].vely = cv[i].ymom / cv[i].dens;
    umin[i].press = dv[i].press;

    umax[i].dens = umin[i].dens;
    umax[i].velx = umin[i].velx;
    umax[i].vely = umin[i].vely;
    umax[i].press = umin[i].press;
  }

  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;

    double rl = cv[i].dens;
    double ul = cv[i].xmom / rl;
    double vl = cv[i].ymom / rl;
    double pl = dv[i].press;

    double rr = cv[j].dens;
    double ur = cv[j].xmom / rr;
    double vr = cv[j].ymom / rr;
    double pr = dv[j].press;

    umin[i].dens = std::min(umin[i].dens, rr);
    umin[i].velx = std::min(umin[i].velx, ur);
    umin[i].vely = std::min(umin[i].vely, vr);
    umin[i].press = std::min(umin[i].press, pr);

    umax[i].dens = std::max(umax[i].dens, rr);
    umax[i].velx = std::max(umax[i].velx, ur);
    umax[i].vely = std::max(umax[i].vely, vr);
    umax[i].press = std::max(umax[i].press, pr);

    umin[j].dens = std::min(umin[j].dens, rl);
    umin[j].velx = std::min(umin[j].velx, ul);
    umin[j].vely = std::min(umin[j].vely, vl);
    umin[j].press = std::min(umin[j].press, pl);

    umax[j].dens = std::max(umax[j].dens, rl);
    umax[j].velx = std::max(umax[j].velx, ul);
    umax[j].vely = std::max(umax[j].vely, vl);
    umax[j].press = std::max(umax[j].press, pl);
  }

  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Periodic &&
        geom.periodicMaster.find(name) != geom.periodicMaster.end()) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        // int i = geom.boundaryNode[ibn].node;
        // int j = geom.boundaryNode[ibn].dummy;
        int i = geom.vertexList[ibn].nodeIdx;
        int j = geom.vertexList[ibn].periodicPair;

        umin[i].dens = std::min(umin[i].dens, umin[j].dens);
        umin[i].velx = std::min(umin[i].velx, umin[j].velx);
        umin[i].vely = std::min(umin[i].vely, umin[j].vely);
        umin[i].press = std::min(umin[i].press, umin[j].press);

        umax[i].dens = std::max(umax[i].dens, umax[j].dens);
        umax[i].velx = std::max(umax[i].velx, umax[j].velx);
        umax[i].vely = std::max(umax[i].vely, umax[j].vely);
        umax[i].press = std::max(umax[i].press, umax[j].press);

        umin[j].dens = umin[i].dens;
        umin[j].velx = umin[i].velx;
        umin[j].vely = umin[i].vely;
        umin[j].press = umin[i].press;

        umax[j].dens = umax[i].dens;
        umax[j].velx = umax[i].velx;
        umax[j].vely = umax[i].vely;
        umax[j].press = umax[i].press;
      }
    }
    ibegn = iendn + 1;
  }
}

void VenkatakrishnanLimiter::limiterUpdate() {
  double eps2[4];
  for (int i = 0; i < geom.phyNodes; ++i) {
    lim[i].dens = 1.0;
    lim[i].velx = 1.0;
    lim[i].vely = 1.0;
    lim[i].press = 1.0;
  }

  double limfac3 = param.limiterCoeff * param.limiterCoeff * param.limiterCoeff;
  double rvolref = 1.0 / std::pow(volRef, 1.5);
  eps2[0] = (limfac3 * rvolref) * (limRef.dens * limRef.dens);
  eps2[1] = (limfac3 * rvolref) * (limRef.velx * limRef.velx);
  eps2[2] = (limfac3 * rvolref) * (limRef.vely * limRef.vely);
  eps2[3] = (limfac3 * rvolref) * (limRef.press * limRef.press);

  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;

    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);
    // Original implementation
    // double voll = std::pow(geom.vol[i], 1.5);
    // double volr = std::pow(geom.vol[j], 1.5);

    // Implementation of Nishikawa in "New Unstructured-Grid Limiter Functions"
    double voll = std::pow(4 * geom.vol[i] / M_PI, 1.5);
    double volr = std::pow(4 * geom.vol[j] / M_PI, 1.5);

    double eps2nl = eps2[0] * voll;
    double d1minl = umin[i].dens - cv[i].dens;
    double d1maxl = umax[i].dens - cv[i].dens;
    double eps2nr = eps2[0] * volr;
    double d1minr = umin[j].dens - cv[j].dens;
    double d1maxr = umax[j].dens - cv[j].dens;
    double d2l = gradx[i].dens * rx + grady[i].dens * ry;
    double d2r = -gradx[j].dens * rx - grady[j].dens * ry;
    double limval = Venkat(d2l, d1minl, d1maxl, eps2nl);
    lim[i].dens = std::min(limval, lim[i].dens);
    limval = Venkat(d2r, d1minr, d1maxr, eps2nr);
    lim[j].dens = std::min(limval, lim[j].dens);

    double ul = cv[i].xmom / cv[i].dens;
    double ur = cv[j].xmom / cv[j].dens;
    eps2nl = eps2[1] * voll;
    d1minl = umin[i].velx - ul;
    d1maxl = umax[i].velx - ul;
    eps2nr = eps2[1] * volr;
    d1minr = umin[j].velx - ur;
    d1maxr = umax[j].velx - ur;
    d2l = gradx[i].velx * rx + grady[i].velx * ry;
    d2r = -gradx[j].velx * rx - grady[j].velx * ry;
    limval = Venkat(d2l, d1minl, d1maxl, eps2nl);
    lim[i].velx = std::min(limval, lim[i].velx);
    limval = Venkat(d2r, d1minr, d1maxr, eps2nr);
    lim[j].velx = std::min(limval, lim[j].velx);

    double vl = cv[i].ymom / cv[i].dens;
    double vr = cv[j].ymom / cv[j].dens;
    eps2nl = eps2[2] * voll;
    d1minl = umin[i].vely - vl;
    d1maxl = umax[i].vely - vl;
    eps2nr = eps2[2] * volr;
    d1minr = umin[j].vely - vr;
    d1maxr = umax[j].vely - vr;
    d2l = gradx[i].vely * rx + grady[i].vely * ry;
    d2r = -gradx[j].vely * rx - grady[j].vely * ry;
    limval = Venkat(d2l, d1minl, d1maxl, eps2nl);
    lim[i].vely = std::min(limval, lim[i].vely);
    limval = Venkat(d2r, d1minr, d1maxr, eps2nr);
    lim[j].vely = std::min(limval, lim[j].vely);

    eps2nl = eps2[3] * voll;
    d1minl = umin[i].press - dv[i].press;
    d1maxl = umax[i].press - dv[i].press;
    eps2nr = eps2[3] * volr;
    d1minr = umin[j].press - dv[j].press;
    d1maxr = umax[j].press - dv[j].press;
    d2l = gradx[i].press * rx + grady[i].press * ry;
    d2r = -gradx[j].press * rx - grady[j].press * ry;
    limval = Venkat(d2l, d1minl, d1maxl, eps2nl);
    lim[i].press = std::min(limval, lim[i].press);
    limval = Venkat(d2r, d1minr, d1maxr, eps2nr);
    lim[j].press = std::min(limval, lim[j].press);
  }

  int ibegn = 0;

  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Periodic &&
        geom.periodicMaster.find(name) != geom.periodicMaster.end()) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        // int i = geom.boundaryNode[ibn].node;
        // int j = geom.boundaryNode[ibn].dummy;
        int i = geom.vertexList[ibn].nodeIdx;
        int j = geom.vertexList[ibn].periodicPair;

        lim[i].dens = std::min(lim[i].dens, lim[j].dens);
        lim[i].velx = std::min(lim[i].velx, lim[j].velx);
        lim[i].vely = std::min(lim[i].vely, lim[j].vely);
        lim[i].press = std::min(lim[i].press, lim[j].press);

        lim[j].dens = lim[i].dens;
        lim[j].velx = lim[i].velx;
        lim[j].vely = lim[i].vely;
        lim[j].press = lim[i].press;
      }
    }

    ibegn = iendn + 1;
  }
}

void NishikawaR3::limiterUpdate() {
  double eps;
  PRIM_VAR init{1.0, 1.0, 1.0, 1.0};
  lim.assign(geom.phyNodes, {1.0, 1.0, 1.0, 1.0});

  double limfac3 = std::pow(param.limiterCoeff, 4);
  double rvolref = 1.0 / std::pow(volRef, 2);
  eps = (limfac3 * rvolref);

  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;

    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);
    // Original implementation
    // double voll = std::pow(geom.vol[i], 1.5);
    // double volr = std::pow(geom.vol[j], 1.5);

    // Implementation of Nishikawa in "New Unstructured-Grid Limiter Functions"
    double voll = std::pow(4 * geom.vol[i] / M_PI, 2);
    double volr = std::pow(4 * geom.vol[j] / M_PI, 2);

    double eps2nl = eps * voll;
    double d1minl = umin[i].dens - cv[i].dens;
    double d1maxl = umax[i].dens - cv[i].dens;
    double eps2nr = eps * volr;
    double d1minr = umin[j].dens - cv[j].dens;
    double d1maxr = umax[j].dens - cv[j].dens;
    double d2l = gradx[i].dens * rx + grady[i].dens * ry;
    double d2r = -gradx[j].dens * rx - grady[j].dens * ry;

    double limval;
    if (d2l >= 0) {
      limval = Nishikawa_R3(d1maxl / limRef.dens, d2l / limRef.dens, eps2nl);
    } else {
      limval = Nishikawa_R3(d1minl / limRef.dens, d2l / limRef.dens, eps2nl);
    }
    lim[i].dens = std::min(limval, lim[i].dens);

    if (d2r >= 0) {
      limval = Nishikawa_R3(d1maxr / limRef.dens, d2r / limRef.dens, eps2nr);
    } else {
      limval = Nishikawa_R3(d1minr / limRef.dens, d2r / limRef.dens, eps2nr);
    }
    lim[j].dens = std::min(limval, lim[j].dens);

    double ul = cv[i].xmom / cv[i].dens;
    double ur = cv[j].xmom / cv[j].dens;
    eps2nl = eps * voll;
    d1minl = umin[i].velx - ul;
    d1maxl = umax[i].velx - ul;
    eps2nr = eps * volr;
    d1minr = umin[j].velx - ur;
    d1maxr = umax[j].velx - ur;
    d2l = gradx[i].velx * rx + grady[i].velx * ry;
    d2r = -gradx[j].velx * rx - grady[j].velx * ry;
    if (d2l >= 0) {
      limval = Nishikawa_R3(d1maxl / limRef.velx, d2l / limRef.velx, eps2nl);
    } else {
      limval = Nishikawa_R3(d1minl / limRef.velx, d2l / limRef.velx, eps2nl);
    }
    lim[i].velx = std::min(limval, lim[i].velx);

    if (d2r >= 0) {
      limval = Nishikawa_R3(d1maxr / limRef.velx, d2r / limRef.velx, eps2nr);
    } else {
      limval = Nishikawa_R3(d1minr / limRef.velx, d2r / limRef.velx, eps2nr);
    }
    lim[j].velx = std::min(limval, lim[j].velx);

    double vl = cv[i].ymom / cv[i].dens;
    double vr = cv[j].ymom / cv[j].dens;
    eps2nl = eps * voll;
    d1minl = umin[i].vely - vl;
    d1maxl = umax[i].vely - vl;
    eps2nr = eps * volr;
    d1minr = umin[j].vely - vr;
    d1maxr = umax[j].vely - vr;
    d2l = gradx[i].vely * rx + grady[i].vely * ry;
    d2r = -gradx[j].vely * rx - grady[j].vely * ry;
    if (d2l >= 0) {
      limval = Nishikawa_R3(d1maxl / limRef.vely, d2l / limRef.vely, eps2nl);
    } else {
      limval = Nishikawa_R3(d1minl / limRef.vely, d2l / limRef.vely, eps2nl);
    }
    lim[i].vely = std::min(limval, lim[i].vely);

    if (d2r >= 0) {
      limval = Nishikawa_R3(d1maxr / limRef.vely, d2r / limRef.vely, eps2nr);
    } else {
      limval = Nishikawa_R3(d1minr / limRef.vely, d2r / limRef.vely, eps2nr);
    }
    lim[j].vely = std::min(limval, lim[j].vely);

    eps2nl = eps * voll;
    d1minl = umin[i].press - dv[i].press;
    d1maxl = umax[i].press - dv[i].press;
    eps2nr = eps * volr;
    d1minr = umin[j].press - dv[j].press;
    d1maxr = umax[j].press - dv[j].press;
    d2l = gradx[i].press * rx + grady[i].press * ry;
    d2r = -gradx[j].press * rx - grady[j].press * ry;
    if (d2l >= 0) {
      limval = Nishikawa_R3(d1maxl / limRef.press, d2l / limRef.press, eps2nl);
    } else {
      limval = Nishikawa_R3(d1minl / limRef.press, d2l / limRef.press, eps2nl);
    }
    lim[i].press = std::min(limval, lim[i].press);

    if (d2r >= 0) {
      limval = Nishikawa_R3(d1maxr / limRef.press, d2r / limRef.press, eps2nr);
    } else {
      limval = Nishikawa_R3(d1minr / limRef.press, d2r / limRef.press, eps2nr);
    }
    lim[j].press = std::min(limval, lim[j].press);
  }

  int ibegn = 0;

  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Periodic &&
        geom.periodicMaster.find(name) != geom.periodicMaster.end()) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        // int i = geom.boundaryNode[ibn].node;
        // int j = geom.boundaryNode[ibn].dummy;
        int i = geom.vertexList[ibn].nodeIdx;
        int j = geom.vertexList[ibn].periodicPair;

        lim[i].dens = std::min(lim[i].dens, lim[j].dens);
        lim[i].velx = std::min(lim[i].velx, lim[j].velx);
        lim[i].vely = std::min(lim[i].vely, lim[j].vely);
        lim[i].press = std::min(lim[i].press, lim[j].press);

        lim[j].dens = lim[i].dens;
        lim[j].velx = lim[i].velx;
        lim[j].vely = lim[i].vely;
        lim[j].press = lim[i].press;
      }
    }

    ibegn = iendn + 1;
  }
}

}  // namespace solver