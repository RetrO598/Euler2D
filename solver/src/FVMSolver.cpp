#include <pre/parameter.h>
#include <solver/FVMSolver.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "solver/limiter.h"
#include "solver/variableDef.h"

namespace solver {

FVMSolver::FVMSolver(preprocess::parameter &parameter,
                     const preprocess::Geometry &geometry)
    : param(parameter), geom(geometry) {
  int nNodes = geom.totNodes;
  cv.reserve(nNodes);
  cvOld.reserve(nNodes);
  diss.reserve(nNodes);
  rhs.reserve(nNodes);

  lim.reserve(nNodes);
  gradx.reserve(nNodes);
  grady.reserve(nNodes);
  umin.reserve(nNodes);
  umax.reserve(nNodes);

  timeSteps.reserve(nNodes);

  if (param.equationtype_ == preprocess::equationType::NavierStokes) {
    gradTx.reserve(nNodes);
    gradTy.reserve(nNodes);
    dvlam.reserve(nNodes);
  } else {
    gradTx.reserve(0);
    gradTy.reserve(0);
    dvlam.reserve(0);
  }

  dv.reserve(nNodes);

  cl = 0.0;

  limiter = new VenkatakrishnanLimiter(param, geom, cv, dv, umin, umax, lim,
                                       gradx, grady);

  rhsIter.reserve(nNodes);
  rhsOld.reserve(nNodes);
  nContr.reserve(nNodes);
}

void FVMSolver::initSolver() {
  double gam1, rgas, temp, mach;

  gam1 = param.gamma - 1.0;
  rgas = gam1 * param.Cp / param.gamma;

  if (param.flowtype_ == preprocess::flowType::External) {
    CONS_VAR init;
    init.dens = param.RhoInfinity;
    init.xmom = param.RhoInfinity * param.uInfinity;
    init.ymom = param.RhoInfinity * param.vInfinity;
    init.ener = param.PsInfinity / gam1 +
                0.5 * param.RhoInfinity * param.velInfinity * param.velInfinity;
    cv.assign(geom.totNodes, init);
  } else {
    double xmin = 1.0e+32;
    double xmax = -1.0e+32;
    for (int i = 0; i < geom.phyNodes; ++i) {
      xmin = std::min(xmin, geom.coords[i].x);
      xmax = std::max(xmax, geom.coords[i].x);
    }
    double dx = xmax - xmin;

    double pinl = param.PsRatio * param.PsOutlet;
    if (pinl >= param.PtInlet) {
      pinl = 0.99999 * param.PtInlet;
    }
    double dp = param.PsOutlet - pinl;
    double dbeta = param.approxFlowAngOut - param.flowAngIn;
    double temp =
        param.TtInlet * std::pow((pinl / param.PtInlet), (gam1 / param.gamma));
    double rho = pinl / (rgas * temp);
    double cs = std::sqrt(param.gamma * pinl / rho);
    double mach = std::sqrt(2.0 * ((param.TtInlet / temp) - 1.0) / gam1);
    double q = mach * cs;

    for (int i = 0; i < geom.totNodes; ++i) {
      double beta = param.flowAngIn + dbeta * (geom.coords[i].x - xmin) / dx;
      double p = pinl + dp * (geom.coords[i].x - xmin) / dx;
      rho = p / (rgas * temp);
      double u = q * std::cos(beta);
      double v = q * std::sin(beta);
      cv[i].dens = rho;
      cv[i].xmom = rho * u;
      cv[i].ymom = rho * v;
      cv[i].ener = p / gam1 + 0.5 * rho * q * q;
    }
  }

  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    if (geom.BoundTypes[ib] >= 700 && geom.BoundTypes[ib] < 800) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.boundNode[ibn].node;
        int j = geom.boundNode[ibn].dummy;
        cv[i].dens = 0.5 * (cv[i].dens + cv[j].dens);
        cv[i].xmom = 0.5 * (cv[i].xmom + cv[j].xmom);
        cv[i].ymom = 0.5 * (cv[i].ymom + cv[j].ymom);
        cv[i].ener = 0.5 * (cv[i].ener + cv[j].ener);
        cv[j].dens = cv[i].dens;
        cv[j].xmom = cv[i].xmom;
        cv[j].ymom = cv[i].ymom;
        cv[j].ener = cv[i].ener;
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::ConvToDependAll() {
  double gam1 = param.gamma - 1.0;
  double rgas = gam1 * param.Cp / param.gamma;
  double g1cp = gam1 * param.Cp;

  if (param.equationtype_ == preprocess::equationType::Euler) {
    for (int i = 0; i < geom.totNodes; ++i) {
      double rhoq = cv[i].xmom * cv[i].xmom + cv[i].ymom * cv[i].ymom;
      dv[i].press = gam1 * (cv[i].ener - 0.5 * rhoq / cv[i].dens);
      dv[i].temp = dv[i].press / (rgas * cv[i].dens);
      dv[i].cs = std::sqrt(g1cp * dv[i].temp);
      dv[i].gamma = param.gamma;
      dv[i].Cp = param.Cp;
    }
  } else {
    double s1 = 110.0;
    double s2 = 288.16;
    double s12 = 1.0 + s1 / s2;
    double cppr = param.Cp / param.Prandtl;
    for (int i = 0; i < geom.totNodes; ++i) {
      double rhoq = cv[i].xmom * cv[i].xmom + cv[i].ymom * cv[i].ymom;
      dv[i].press = gam1 * (cv[i].ener - 0.5 * rhoq / cv[i].dens);
      dv[i].temp = dv[i].press / (rgas * cv[i].dens);
      dv[i].cs = std::sqrt(g1cp * dv[i].temp);
      dv[i].gamma = param.gamma;
      dv[i].Cp = param.Cp;
      double rat = std::sqrt(dv[i].temp / s2) * s12 / (1.0 + s1 / dv[i].temp);
      dvlam[i].mu = param.refVisc * rat;
      dvlam[i].lambda = dvlam[i].mu * cppr;
    }
  }
}

void FVMSolver::ConvToDepend(int i) {
  double gam1 = param.gamma - 1.0;
  double rgas = gam1 * param.Cp / param.gamma;
  double g1cp = gam1 * param.Cp;

  if (param.equationtype_ == preprocess::equationType::Euler) {
    double rhoq = cv[i].xmom * cv[i].xmom + cv[i].ymom * cv[i].ymom;
    dv[i].press = gam1 * (cv[i].ener - 0.5 * rhoq / cv[i].dens);
    dv[i].temp = dv[i].press / (rgas * cv[i].dens);
    dv[i].cs = std::sqrt(g1cp * dv[i].temp);
    dv[i].gamma = param.gamma;
    dv[i].Cp = param.Cp;
  } else {
    double s1 = 110.0;
    double s2 = 288.16;
    double s12 = 1.0 + s1 / s2;

    double rhoq = cv[i].xmom * cv[i].xmom + cv[i].ymom * cv[i].ymom;
    dv[i].press = gam1 * (cv[i].ener - 0.5 * rhoq / cv[i].dens);
    dv[i].temp = dv[i].press / (rgas * cv[i].dens);
    dv[i].cs = std::sqrt(g1cp * dv[i].temp);
    dv[i].gamma = param.gamma;
    dv[i].Cp = param.Cp;
    double rat = std::sqrt(dv[i].temp / s2) * s12 / (1.0 + s1 / dv[i].temp);
    dvlam[i].mu = param.refVisc * rat;
    dvlam[i].lambda = dvlam[i].mu * (param.Cp / param.Prandtl);
  }
}
void FVMSolver::DissipInit(int irk, double beta) {
  if (irk == 0 || beta > 0.99) {
    for (int i = 0; i < geom.totNodes; ++i) {
      diss[i].dens = 0.0;
      diss[i].xmom = 0.0;
      diss[i].ymom = 0.0;
      diss[i].ener = 0.0;
    }
  } else {
    double blend = 1.0 - beta;
    for (int i = 0; i < geom.totNodes; ++i) {
      diss[i].dens *= blend;
      diss[i].xmom *= blend;
      diss[i].ymom *= blend;
      diss[i].ener *= blend;
    }
  }
}

void FVMSolver::DissipRoe2(double beta) {
  double beta5 = 0.5 * beta;
  double rrho, gam1, ggm1, rl, ul, vl, pl, hl, rr, ur, vr, pr, hr, rav, dd, dd1,
      uav, vav, hav, q2a, c2a, cav, uv, du, h1, h2, h3, h4, h5, delta, eabs1,
      eabs2, eabs4;
  double fd[4];

  for (int ie = 0; ie < geom.totEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double nx = geom.sij[ie].x / ds;
    double ny = geom.sij[ie].y / ds;

    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);

    rrho = 1.0 / cv[i].dens;
    gam1 = dv[i].gamma - 1.0;
    ggm1 = dv[i].gamma / gam1;
    rl = cv[i].dens + lim[i].dens * (gradx[i].dens * rx + grady[i].dens * ry);
    ul = cv[i].xmom * rrho +
         lim[i].velx * (gradx[i].velx * rx + grady[i].velx * ry);
    vl = cv[i].ymom * rrho +
         lim[i].vely * (gradx[i].vely * rx + grady[i].vely * ry);
    pl = dv[i].press +
         lim[i].press * (gradx[i].press * rx + grady[i].press * ry);
    hl = ggm1 * pl / rl + 0.5 * (ul * ul + vl * vl);

    rrho = 1.0 / cv[j].dens;
    gam1 = dv[j].gamma - 1.0;
    ggm1 = dv[j].gamma / gam1;
    rr = cv[j].dens - lim[j].dens * (gradx[j].dens * rx + grady[j].dens * ry);
    ur = cv[j].xmom * rrho -
         lim[j].velx * (gradx[j].velx * rx + grady[j].velx * ry);
    vr = cv[j].ymom * rrho -
         lim[j].vely * (gradx[j].vely * rx + grady[j].vely * ry);
    pr = dv[j].press -
         lim[j].press * (gradx[j].press * rx + grady[j].press * ry);
    hr = ggm1 * pr / rr + 0.5 * (ur * ur + vr * vr);

    rav = std::sqrt(rl * rr);
    gam1 = 0.5 * (dv[i].gamma + dv[j].gamma) - 1.0;
    dd = rav / rl;
    dd1 = 1.0 / (1.0 + dd);
    uav = (ul + dd * ur) * dd1;
    vav = (vl + dd * vr) * dd1;
    hav = (hl + dd * hr) * dd1;
    q2a = 0.5 * (uav * uav + vav * vav);
    c2a = gam1 * (hav - q2a);
    cav = std::sqrt(c2a);
    uv = uav * nx + vav * ny;
    du = (ur - ul) * nx + (vr - vl) * ny;

    h1 = std::abs(uv - cav);
    h2 = std::abs(uv);
    h4 = std::abs(uv + cav);
    delta = param.entropyCorreRoe * h4;

    eabs1 = EntropyCorr(h1, delta);
    eabs2 = EntropyCorr(h2, delta);
    eabs4 = EntropyCorr(h4, delta);

    h1 = rav * cav * du;
    h2 = eabs1 * (pr - pl - h1) / (2.0 * c2a);
    h3 = eabs2 * (rr - rl - (pr - pl) / c2a);
    h4 = eabs2 * rav;
    h5 = eabs4 * (pr - pl + h1) / (2.0 * c2a);

    fd[0] = h2 + h3 + h5;
    fd[1] = h2 * (uav - cav * nx) + h3 * uav + h4 * (ur - ul - du * nx) +
            h5 * (uav + cav * nx);
    fd[2] = h2 * (vav - cav * ny) + h3 * vav + h4 * (vr - vl - du * ny) +
            h5 * (vav + cav * ny);
    fd[3] = h2 * (hav - cav * uv) + h3 * q2a +
            h4 * (uav * (ur - ul) + vav * (vr - vl) - uv * du) +
            h5 * (hav + cav * uv);

    ds *= beta5;
    diss[i].dens += fd[0] * ds;
    diss[i].xmom += fd[1] * ds;
    diss[i].ymom += fd[2] * ds;
    diss[i].ener += fd[3] * ds;

    diss[j].dens -= fd[0] * ds;
    diss[j].xmom -= fd[1] * ds;
    diss[j].ymom -= fd[2] * ds;
    diss[j].ener -= fd[3] * ds;
  }
}

double FVMSolver::EntropyCorr(double z, double d) {
  if (z > d) {
    return z;
  } else {
    return 0.5 * (z * z + d * d) / d;
  }
}

void FVMSolver::FluxRoe2() {
  double rx, ry, rrho, gam1, ggm1, rl, ul, vl, pl, hl, rr, ur, vr, pr, hr, qsrl,
      qsrr, pav;
  double fc[4];

  for (int i = 0; i < geom.totNodes; ++i) {
    rhs[i].dens = -diss[i].dens;
    rhs[i].xmom = -diss[i].xmom;
    rhs[i].ymom = -diss[i].ymom;
    rhs[i].ener = -diss[i].ener;
  }

  for (int ie = 0; ie < geom.totEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    ry = 0.5 * (geom.coords[i].y - geom.coords[i].y);

    rrho = 1.0 / cv[i].dens;
    gam1 = dv[i].gamma - 1.0;
    ggm1 = dv[i].gamma / gam1;
    rl = cv[i].dens + lim[i].dens * (gradx[i].dens * rx + grady[i].dens * ry);
    ul = cv[i].xmom * rrho +
         lim[i].velx * (gradx[i].velx * rx + grady[i].velx * ry);
    vl = cv[i].ymom * rrho +
         lim[i].vely * (gradx[i].vely * rx + grady[i].vely * ry);
    pl = dv[i].press +
         lim[i].press * (gradx[i].press * rx + grady[i].press * ry);
    hl = ggm1 * pl / rl + 0.5 * (ul * ul + vl * vl);
    qsrl = (ul * geom.sij[ie].x + vl * geom.sij[ie].y) * rl;

    rrho = 1.0 / cv[j].dens;
    gam1 = dv[j].gamma - 1.0;
    ggm1 = dv[j].gamma / gam1;
    rr = cv[j].dens - lim[j].dens * (gradx[j].dens * rx + grady[j].dens * ry);
    ur = cv[j].xmom * rrho -
         lim[j].velx * (gradx[j].velx * rx + grady[j].velx * ry);
    vr = cv[j].ymom * rrho -
         lim[j].vely * (gradx[j].vely * rx + grady[j].vely * ry);
    pr = dv[j].press -
         lim[j].press * (gradx[j].press * rx + grady[j].press * ry);
    hr = ggm1 * pr / rr + 0.5 * (ur * ur + vr * vr);
    qsrr = (ur * geom.sij[ie].x + vr * geom.sij[ie].y) * rr;

    pav = 0.5 * (pl + pr);
    fc[0] = 0.5 * (qsrl + qsrr);
    fc[1] = 0.5 * (qsrl * ul + qsrr * ur) + geom.sij[ie].x * pav;
    fc[2] = 0.5 * (qsrl * vl + qsrr * vr) + geom.sij[ie].y * pav;
    fc[3] = 0.5 * (qsrl * hl + qsrr * hr);

    rhs[i].dens += fc[0];
    rhs[i].xmom += fc[1];
    rhs[i].ymom += fc[2];
    rhs[i].ener += fc[3];

    rhs[j].dens -= fc[0];
    rhs[j].xmom -= fc[1];
    rhs[j].ymom -= fc[2];
    rhs[j].ener -= fc[3];
  }

  FluxWalls();
}

void FVMSolver::FluxWalls() {
  int ibegf = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendf = geom.ibound[ib].bfaceIndex;
    if (geom.BoundTypes[ib] >= 300 && geom.BoundTypes[ib] < 500) {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        int i = geom.boundFace[ibf].nodei;
        int j = geom.boundFace[ibf].nodej;
        double sx = geom.sbf[ibf].x / 12.0;
        double sy = geom.sbf[ibf].y / 12.0;
        double pl = 5.0 * dv[i].press + dv[j].press;
        double pr = dv[i].press + 5.0 * dv[j].press;
        rhs[i].xmom += sx * pl;
        rhs[i].ymom += sy * pl;
        rhs[j].xmom += sx * pr;
        rhs[j].ymom += sy * pr;
      }
    }
    ibegf = iendf + 1;
  }
}
}  // namespace solver