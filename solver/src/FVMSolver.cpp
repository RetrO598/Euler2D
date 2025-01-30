#include "solver/numeric.h"
#include <pre/parameter.h>
#include <solver/FVMSolver.h>
#include <solver/limiter.h>
#include <solver/variableDef.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

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

  limiter = new NishikawaR3(param, geom, cv, dv, umin, umax, lim, gradx, grady);

  numeric = new NumericSLAU2(param, geom, cv, dv, diss, rhs, lim, gradx, grady);

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
        int i = geom.boundaryNode[ibn].node;
        int j = geom.boundaryNode[ibn].dummy;
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
} // namespace solver