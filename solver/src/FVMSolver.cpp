#include <pre/parameter.h>
#include <solver/FVMSolver.h>
#include <solver/limiter.h>
#include <solver/numeric.h>
#include <solver/variableDef.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <vector>

#include "solver/timeIntegrator.hpp"

namespace solver {

FVMSolver::FVMSolver(preprocess::parameter &parameter,
                     const preprocess::Geometry &geometry)
    : param(parameter), geom(geometry) {
  int nNodes = geom.phyNodes;
  cv.resize(nNodes);
  cvOld.resize(nNodes);
  diss.resize(nNodes);
  rhs.resize(nNodes);

  lim.resize(nNodes);
  gradx.resize(nNodes);
  grady.resize(nNodes);
  umin.resize(nNodes);
  umax.resize(nNodes);

  timeSteps.resize(nNodes);

  if (param.equationtype_ == preprocess::equationType::NavierStokes) {
    gradTx.resize(nNodes);
    gradTy.resize(nNodes);
    dvlam.resize(nNodes);
  } else {
    gradTx.resize(0);
    gradTy.resize(0);
    dvlam.resize(0);
  }

  diag.resize(nNodes);
  intermediateSol.resize(nNodes);
  increment.resize(nNodes);

  dv.resize(nNodes);

  cl = 0.0;

  if (param.convecScheme == preprocess::ConvectionScheme::ROE) {
    std::cout << "using ROE scheme." << "\n";
    numeric = std::make_unique<NumericRoe>(param, geom, cv, dv, diss, rhs, lim,
                                           gradx, grady);
  } else if (param.convecScheme == preprocess::ConvectionScheme::SLAU2) {
    std::cout << "using SLAU2 scheme." << "\n";
    numeric = std::make_unique<NumericSLAU2>(param, geom, cv, dv, diss, rhs,
                                             lim, gradx, grady);
  } else if (param.convecScheme == preprocess::ConvectionScheme::AUSM) {
    std::cout << "using AUSM scheme." << "\n";
    numeric = std::make_unique<NumericAUSM>(param, geom, cv, dv, diss, rhs, lim,
                                            gradx, grady);
  } else if (param.convecScheme == preprocess::ConvectionScheme::AUSMUP2) {
    std::cout << "using AUSMUP2 scheme." << "\n";
    numeric = std::make_unique<NumericAUSMUP2>(param, geom, cv, dv, diss, rhs,
                                               lim, gradx, grady);
  }

  if (param.limiterType == preprocess::Limiter::VenkataKrishnan) {
    std::cout << "using VenkataKrishnan limiter." << "\n";
    limiter = std::make_unique<VenkatakrishnanLimiter>(
        param, geom, cv, dv, umin, umax, lim, gradx, grady);
  } else if (param.limiterType == preprocess::Limiter::NishikawaR3) {
    std::cout << "using NishikawaR3 limiter." << "\n";
    limiter = std::make_unique<NishikawaR3>(param, geom, cv, dv, umin, umax,
                                            lim, gradx, grady);
  }

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
    cv.assign(geom.phyNodes, init);
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
    for (int i = 0; i < geom.phyNodes; ++i) {
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
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Periodic) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.vertexList[ibn].nodeIdx;
        int j = geom.vertexList[ibn].periodicPair;
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
    for (int i = 0; i < geom.phyNodes; ++i) {
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
    for (int i = 0; i < geom.phyNodes; ++i) {
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

void FVMSolver::updateResidualRK(int irk) {
  numeric->DissipInit();

  if (param.equationtype_ == preprocess::equationType::NavierStokes) {
    GradientsVisc();
    fluxVisc();
  } else if (param.equationtype_ == preprocess::equationType::Euler) {
    Gradients();
  }

  limiter->limiterInit();
  limiter->limiterUpdate();

  numeric->FluxNumeric();
  BoundaryConditions();
  ZeroRes();
  PeriodicCons(rhs);

  double fac = param.stageCoeff[irk] * param.CFL;
  for (int i = 0; i < geom.phyNodes; ++i) {
    double adtv = fac * timeSteps[i] / geom.vol[i];
    rhs[i].dens *= adtv;
    rhs[i].xmom *= adtv;
    rhs[i].ymom *= adtv;
    rhs[i].ener *= adtv;
  }
}

void FVMSolver::assignCVold() { cvOld = cv; }

void FVMSolver::updateCV() {
  for (int i = 0; i < geom.phyNodes; ++i) {
    cv[i].dens = cvOld[i].dens - rhs[i].dens;
    cv[i].xmom = cvOld[i].xmom - rhs[i].xmom;
    cv[i].ymom = cvOld[i].ymom - rhs[i].ymom;
    cv[i].ener = cvOld[i].ener - rhs[i].ener;
  }
}

void FVMSolver::computeWaveSpeed() {
  std::size_t sum = 0;
  for (auto &i : geom.pointList) {
    sum += i.getnPoint();
  }

  waveSpeed.resize(sum);
  waveSpeedJ.resize(sum);

  std::size_t index = 0;
  for (std::size_t iPoint = 0; iPoint < geom.pointList.size(); ++iPoint) {
    for (std::size_t j = 0; j < geom.pointList[iPoint].getnPoint(); ++j) {
      auto jPoint = geom.pointList[iPoint].getPoint(j);
      auto jEdge = geom.pointList[iPoint].getEdge(j);

      double sx = geom.sij[jEdge].x;
      double sy = geom.sij[jEdge].y;

      double ds = std::sqrt(sx * sx + sy * sy);

      double nx = sx / ds;
      double ny = sy / ds;

      double velx = 0.5 * (cv[iPoint].xmom / cv[iPoint].dens +
                           cv[jPoint].xmom / cv[jPoint].dens);
      double vely = 0.5 * (cv[iPoint].ymom / cv[iPoint].dens +
                           cv[jPoint].ymom / cv[jPoint].dens);

      double cs = 0.5 * (dv[iPoint].cs + dv[jPoint].cs);

      double wavespeed = (std::abs(velx * nx + vely * ny) + cs) * ds;

      double uj = cv[jPoint].xmom / cv[jPoint].dens;
      double vj = cv[jPoint].ymom / cv[jPoint].dens;

      double wavespeedj = (std::abs(uj * nx + vj * ny) + dv[jPoint].cs) * ds;

      waveSpeedJ[index] = wavespeedj;
      waveSpeed[index] = wavespeed;
      index++;
    }
  }
}

void FVMSolver::computeJacobianDiag(double factor) {
  std::size_t index = 0;
  for (std::size_t iPoint = 0; iPoint < geom.pointList.size(); ++iPoint) {
    diag[iPoint] = geom.vol[iPoint] / (timeSteps[iPoint] * param.CFL);
    for (std::size_t j = 0; j < geom.pointList[iPoint].getnPoint(); ++j) {
      double lambda = waveSpeed[index] * factor;
      waveSpeed[index] = lambda;
      diag[iPoint] += 0.5 * lambda;
      index++;
    }
  }
}

void FVMSolver::lowerSweep() {
  CONS_VAR df{0.0, 0.0, 0.0, 0.0};

  std::size_t index = 0;
  for (std::size_t iPoint = 0; iPoint < geom.pointList.size(); ++iPoint) {
    df = {0.0, 0.0, 0.0, 0.0};
    for (std::size_t ineighbor = 0;
         ineighbor < geom.pointList[iPoint].getnPoint(); ++ineighbor) {
      std::size_t jPoint = geom.pointList[iPoint].getPoint(ineighbor);
      std::size_t jEdge = geom.pointList[iPoint].getEdge(ineighbor);

      if (jPoint < iPoint) {
        CONS_VAR sol = cv[jPoint];
        CONS_VAR dSol;
        dSol.dens = sol.dens + intermediateSol[jPoint].dens;
        dSol.xmom = sol.xmom + intermediateSol[jPoint].xmom;
        dSol.ymom = sol.ymom + intermediateSol[jPoint].ymom;
        dSol.ener = sol.ener + intermediateSol[jPoint].ener;

        double sx, sy, nx, ny, gamma;

        sx = geom.sij[jEdge].x;
        sy = geom.sij[jEdge].y;
        double ds = std::sqrt(sx * sx + sy * sy);
        nx = sx / ds;
        ny = sy / ds;

        if (geom.edge[jEdge].nodei == jPoint) {
          nx = -nx;
          ny = -ny;
        }

        gamma = dv[jPoint].gamma;

        auto f = fluxDifference::computeDifference(sol, nx, ny, gamma);
        auto dfj = fluxDifference::computeDifference(dSol, nx, ny, gamma);

        dfj.dens -= f.dens;
        dfj.xmom -= f.xmom;
        dfj.ymom -= f.ymom;
        dfj.ener -= f.ener;

        df.dens +=
            (dfj.dens * ds - waveSpeed[index] * intermediateSol[jPoint].dens);
        df.xmom +=
            (dfj.xmom * ds - waveSpeed[index] * intermediateSol[jPoint].xmom);
        df.ymom +=
            (dfj.ymom * ds - waveSpeed[index] * intermediateSol[jPoint].ymom);
        df.ener +=
            (dfj.ener * ds - waveSpeed[index] * intermediateSol[jPoint].ener);
      }
      index++;
    }

    intermediateSol[iPoint].dens =
        (-rhs[iPoint].dens - 0.5 * df.dens) / diag[iPoint];
    intermediateSol[iPoint].xmom =
        (-rhs[iPoint].xmom - 0.5 * df.xmom) / diag[iPoint];
    intermediateSol[iPoint].ymom =
        (-rhs[iPoint].ymom - 0.5 * df.ymom) / diag[iPoint];
    intermediateSol[iPoint].ener =
        (-rhs[iPoint].ener - 0.5 * df.ener) / diag[iPoint];
  }
}

void FVMSolver::upperSweep() {
  CONS_VAR df{0.0, 0.0, 0.0, 0.0};

  std::size_t index = waveSpeed.size();
  for (int iPoint = geom.pointList.size() - 1; iPoint > -1; --iPoint) {
    df = {0.0, 0.0, 0.0, 0.0};
    for (int ineighbor = geom.pointList[iPoint].getnPoint() - 1; ineighbor > -1;
         --ineighbor) {
      index--;
      std::size_t jPoint = geom.pointList[iPoint].getPoint(ineighbor);
      std::size_t jEdge = geom.pointList[iPoint].getEdge(ineighbor);

      if (jPoint > iPoint) {
        CONS_VAR sol = cv[jPoint];
        CONS_VAR dSol;
        dSol.dens = sol.dens + increment[jPoint].dens;
        dSol.xmom = sol.xmom + increment[jPoint].xmom;
        dSol.ymom = sol.ymom + increment[jPoint].ymom;
        dSol.ener = sol.ener + increment[jPoint].ener;

        double sx, sy, nx, ny, gamma;

        sx = geom.sij[jEdge].x;
        sy = geom.sij[jEdge].y;
        double ds = std::sqrt(sx * sx + sy * sy);
        nx = sx / ds;
        ny = sy / ds;

        if (geom.edge[jEdge].nodei == jPoint) {
          nx = -nx;
          ny = -ny;
        }

        gamma = dv[jPoint].gamma;

        auto f = fluxDifference::computeDifference(sol, nx, ny, gamma);
        auto dfj = fluxDifference::computeDifference(dSol, nx, ny, gamma);

        dfj.dens -= f.dens;
        dfj.xmom -= f.xmom;
        dfj.ymom -= f.ymom;
        dfj.ener -= f.ener;

        df.dens += (dfj.dens * ds - waveSpeed[index] * increment[jPoint].dens);
        df.xmom += (dfj.xmom * ds - waveSpeed[index] * increment[jPoint].xmom);
        df.ymom += (dfj.ymom * ds - waveSpeed[index] * increment[jPoint].ymom);
        df.ener += (dfj.ener * ds - waveSpeed[index] * increment[jPoint].ener);
      }
    }

    increment[iPoint].dens =
        intermediateSol[iPoint].dens - 0.5 * df.dens / diag[iPoint];
    increment[iPoint].xmom =
        intermediateSol[iPoint].xmom - 0.5 * df.xmom / diag[iPoint];
    increment[iPoint].ymom =
        intermediateSol[iPoint].ymom - 0.5 * df.ymom / diag[iPoint];
    increment[iPoint].ener =
        intermediateSol[iPoint].ener - 0.5 * df.ener / diag[iPoint];
  }
}

void FVMSolver::LUSGSupdate() {
  for (std::size_t i = 0; i < geom.phyNodes; ++i) {
    cv[i].dens += increment[i].dens;
    cv[i].xmom += increment[i].xmom;
    cv[i].ymom += increment[i].ymom;
    cv[i].ener += increment[i].ener;
  }
}

void FVMSolver::computeResidualLUSGS() {
  numeric->DissipInit();

  if (param.equationtype_ == preprocess::equationType::NavierStokes) {
    GradientsVisc();
    fluxVisc();
  } else if (param.equationtype_ == preprocess::equationType::Euler) {
    Gradients();
  }

  limiter->limiterInit();
  limiter->limiterUpdate();

  numeric->FluxNumeric();
  BoundaryConditions();
  ZeroRes();
  PeriodicCons(rhs);
}
}  // namespace solver