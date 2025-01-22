#include <pre/parameter.h>
#include <solver/FVMSolver.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
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

  limiter = std::make_unique<VenkatakrishnanLimiter>();

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

  // std::ofstream output("init.txt");
  // for (int i = 0; i < geom.totNodes; ++i) {
  //   output << cv[i].dens << " " << cv[i].xmom << " " << cv[i].ymom << " "
  //          << cv[i].ener << "\n";
  // }
  // // std::cout << "init end" << "\n";
  // output.close();
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

  // std::ofstream output("depend.txt");
  // for (int i = 0; i < geom.totNodes; ++i) {
  //   output << dv[i].press << " " << dv[i].temp << " " << dv[i].cs << " "
  //          << dv[i].gamma << " " << dv[i].Cp << "\n";
  // }
  // // std::cout << "init end" << "\n";
  // output.close();
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

void FVMSolver::BoundaryConditions() {
  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int itype = geom.BoundTypes[ib];
    int iendn = geom.ibound[ib].bnodeIndex;

    if (itype >= 100 && itype < 200) {
      BoundInflow(ibegn, iendn);
    } else if (itype >= 200 && itype < 300) {
      BoundOutflow(ibegn, iendn);
    } else if (itype >= 600 && itype < 700) {
      BoundFarfield(ibegn, iendn);
    }
    ibegn = iendn + 1;
  }

  ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int itype = geom.BoundTypes[ib];
    int iendn = geom.ibound[ib].bnodeIndex;
    if (param.equationtype_ == preprocess::equationType::NavierStokes &&
        (itype >= 300 && itype < 400)) {
      BoundWallVisc(ibegn, iendn);
    }
    ibegn = iendn + 1;
  }

  std::ofstream outputCV("cv.txt");
  for (int i = 0; i < geom.totNodes; ++i) {
    outputCV << cv[i].dens << " " << cv[i].xmom << " " << cv[i].ymom << " "
             << cv[i].ener << "\n";
  }
  // std::cout << "init end" << "\n";
  outputCV.close();

  std::ofstream outputDV("depend.txt");
  for (int i = 0; i < geom.totNodes; ++i) {
    outputDV << dv[i].press << " " << dv[i].temp << " " << dv[i].cs << " "
             << dv[i].gamma << " " << dv[i].Cp << "\n";
  }
  // std::cout << "init end" << "\n";
  outputDV.close();
}

void FVMSolver::BoundInflow(int beg, int end) {
  for (int ib = beg; ib <= end; ++ib) {
    int ibn = geom.boundNode[ib].node;
    int idn = geom.boundNode[ib].dummy;
    int ie = geom.boundNode[ib].indexEdge;

    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double sxn = geom.sij[ie].x / ds;
    double syn = geom.sij[ie].y / ds;

    double gam1 = dv[ibn].gamma - 1.0;
    double ggm1 = dv[ibn].gamma / gam1;
    double rgas = gam1 * dv[ibn].Cp / dv[ibn].gamma;
    double rrho = 1.0 / cv[ibn].dens;
    double u = cv[ibn].xmom * rrho;
    double v = cv[ibn].ymom * rrho;
    double uabs = std::sqrt(u * u + v * v);
    double unorm = u * sxn + v * syn;
    double cosa = 0.0;
    if (uabs < 1e-16) {
      cosa = 1.0;
    } else {
      cosa = -unorm / uabs;
    }

    double c02 = dv[ibn].cs * dv[ibn].cs + 0.5 * gam1 * uabs * uabs;
    double rinv = unorm - 2.0 * dv[ibn].cs / gam1;
    double dis =
        (gam1 * cosa * cosa + 2.0) * c02 / (gam1 * rinv * rinv) - 0.5 * gam1;
    dis = std::max(dis, 1e-16);
    double cb = -rinv * (gam1 / (gam1 * cosa * cosa + 2.0)) *
                (1.0 + cosa * std::sqrt(dis));
    double cc02 = std::min((cb * cb) / c02, 1.0);
    double tb = param.TtInlet * cc02;
    double pb = param.PtInlet * std::pow((tb / param.TtInlet), ggm1);
    double rhob = pb / (rgas * tb);
    double uabsb = (2.0 * dv[ibn].Cp) * (param.TtInlet - tb);
    uabsb = std::sqrt(uabsb);
    double ub = uabsb * std::cos(param.flowAngIn);
    double vb = uabsb * std::sin(param.flowAngIn);

    cv[idn].dens = rhob;
    cv[idn].xmom = rhob * ub;
    cv[idn].ymom = rhob * vb;
    cv[idn].ener = pb / gam1 + 0.5 * rhob * (ub * ub + vb * vb);

    ConvToDepend(idn);
  }
}

void FVMSolver::BoundFarfield(int beg, int end) {
  PRIM_VAR farfield;
  if (param.farCorrect) {
    double bet = std::sqrt(1.0 - param.MaInfinity * param.MaInfinity);
    double cir = 0.25 * param.chord * cl * param.velInfinity / M_PI;
    for (int ib = beg; ib <= end; ++ib) {
      int ibn = geom.boundNode[ib].node;
      int idn = geom.boundNode[ib].dummy;
      int ie = geom.boundNode[ib].indexEdge;

      double gam1 = dv[ibn].gamma - 1.0;
      double ggm1 = dv[ibn].gamma / gam1;
      double gmr = 1.0 / dv[ibn].gamma;
      double gmg = gam1 / dv[ibn].gamma;
      double xa = geom.coords[ibn].x - param.xRefPoint;
      double ya = geom.coords[ibn].y - param.yRefPoint;
      double dist = std::sqrt(xa * xa + ya * ya);
      double angle = std::atan2(ya, xa);
      double sn = std::sin(angle - param.Aoa);
      double dn = 1.0 - param.MaInfinity * param.MaInfinity * sn * sn;
      double vc = cir * bet / (dn * dist);

      farfield.velx = param.uInfinity + vc * std::sin(angle);
      farfield.vely = param.vInfinity - vc * std::cos(angle);
      double qv2 =
          farfield.velx * farfield.velx + farfield.vely * farfield.vely;
      farfield.press =
          std::pow(std::pow(param.PsInfinity, gmg) +
                       gmg * param.RhoInfinity *
                           (param.velInfinity * param.velInfinity - qv2) /
                           (2.0 * std::pow(param.PsInfinity, gmr)),
                   ggm1);
      farfield.dens =
          param.RhoInfinity * std::pow(farfield.press / param.PsInfinity, gmr);

      double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                            geom.sij[ie].y * geom.sij[ie].y);
      double sxn = geom.sij[ie].x / ds;
      double syn = geom.sij[ie].y / ds;

      double rhoe = cv[ibn].dens;
      double ue = cv[ibn].xmom / rhoe;
      double ve = cv[ibn].ymom / rhoe;
      double pe = dv[ibn].press;
      double qn = sxn * ue + syn * ve;
      double crho0 = dv[ibn].cs * rhoe;
      double rhoa;
      double ua;
      double va;
      double pa;
      double sgn;
      double pb;

      if (param.MaInfinity < 1.0) {
        if (qn < 0.0) {
          rhoa = farfield.dens;
          ua = farfield.velx;
          va = farfield.vely;
          pa = farfield.press;
          sgn = -1.0;
          pb = 0.5 * (pa + pe - crho0 * (sxn * (ua - ue) + syn * (va - ve)));
        } else {
          rhoa = rhoe;
          ua = ue;
          va = ve;
          pa = pe;
          sgn = 1.0;
          pb = farfield.press;
        }
        double dprhoc = sgn * (pa - pb) / crho0;
        cv[idn].dens = rhoa + (pb - pa) / (dv[ibn].cs * dv[ibn].cs);
        cv[idn].xmom = cv[idn].dens * (ua + sxn * dprhoc);
        cv[idn].ymom = cv[idn].dens * (va + syn * dprhoc);
        cv[idn].ener =
            pb / gam1 +
            0.5 * (cv[idn].xmom * cv[idn].xmom + cv[idn].ymom * cv[idn].ymom) /
                cv[idn].dens;
      } else {
        if (qn < 0.0) {
          cv[idn].dens = param.RhoInfinity;
          cv[idn].xmom = param.RhoInfinity * param.uInfinity;
          cv[idn].ymom = param.RhoInfinity * param.vInfinity;
          cv[idn].ener = param.PsInfinity / gam1 + 0.5 * param.RhoInfinity *
                                                       param.velInfinity *
                                                       param.velInfinity;
        } else {
          cv[idn].dens = cv[ibn].dens;
          cv[idn].xmom = cv[ibn].xmom;
          cv[idn].ymom = cv[ibn].ymom;
          cv[idn].ener = cv[ibn].ener;
        }
      }
      ConvToDepend(idn);
    }
  } else {
    for (int ib = beg; ib <= end; ++ib) {
      int ibn = geom.boundNode[ib].node;
      int idn = geom.boundNode[ib].dummy;
      int ie = geom.boundNode[ib].indexEdge;

      double gam1 = dv[ibn].gamma - 1.0;
      farfield.dens = param.RhoInfinity;
      farfield.velx = param.uInfinity;
      farfield.vely = param.vInfinity;
      farfield.press = param.PsInfinity;

      double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                            geom.sij[ie].y * geom.sij[ie].y);
      double sxn = geom.sij[ie].x / ds;
      double syn = geom.sij[ie].y / ds;

      double rhoe = cv[ibn].dens;
      double ue = cv[ibn].xmom / rhoe;
      double ve = cv[ibn].ymom / rhoe;
      double pe = dv[ibn].press;
      double qn = sxn * ue + syn * ve;
      double crho0 = dv[ibn].cs * rhoe;
      double rhoa;
      double ua;
      double va;
      double pa;
      double sgn;
      double pb;

      if (param.MaInfinity < 1.0) {
        if (qn < 0.0) {
          rhoa = farfield.dens;
          ua = farfield.velx;
          va = farfield.vely;
          pa = farfield.press;
          sgn = -1.0;
          pb = 0.5 * (pa + pe - crho0 * (sxn * (ua - ue) + syn * (va - ve)));
        } else {
          rhoa = rhoe;
          ua = ue;
          va = ve;
          pa = pe;
          sgn = 1.0;
          pb = farfield.press;
        }
        double dprhoc = sgn * (pa - pb) / crho0;
        cv[idn].dens = rhoa + (pb - pa) / (dv[ibn].cs * dv[ibn].cs);
        cv[idn].xmom = cv[idn].dens * (ua + sxn * dprhoc);
        cv[idn].ymom = cv[idn].dens * (va + syn * dprhoc);
        cv[idn].ener =
            pb / gam1 +
            0.5 * (cv[idn].xmom * cv[idn].xmom + cv[idn].ymom * cv[idn].ymom) /
                cv[idn].dens;
      } else {
        if (qn < 0.0) {
          cv[idn].dens = param.RhoInfinity;
          cv[idn].xmom = param.RhoInfinity * param.uInfinity;
          cv[idn].ymom = param.RhoInfinity * param.vInfinity;
          cv[idn].ener = param.PsInfinity / gam1 + 0.5 * param.RhoInfinity *
                                                       param.velInfinity *
                                                       param.velInfinity;
        } else {
          cv[idn].dens = cv[ibn].dens;
          cv[idn].xmom = cv[ibn].xmom;
          cv[idn].ymom = cv[ibn].ymom;
          cv[idn].ener = cv[ibn].ener;
        }
      }
      ConvToDepend(idn);
    }
  }
}

void FVMSolver::BoundOutflow(int beg, int end) {
  for (int ib = beg; ib <= end; ++ib) {
    int ibn = geom.boundNode[ib].node;
    int idn = geom.boundNode[ib].dummy;
    int ie = geom.boundNode[ib].indexEdge;

    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double sxn = geom.sij[ie].x / ds;
    double syn = geom.sij[ie].y / ds;

    double gam1 = dv[ibn].gamma - 1.0;
    double rrho = 1.0 / cv[ibn].dens;
    double u = cv[ibn].xmom * rrho;
    double v = cv[ibn].ymom * rrho;
    double q = std::sqrt(u * u + v * v);
    double mach = q / dv[ibn].cs;

    if (mach < 1.0) {
      double rrhoc = rrho / dv[ibn].cs;
      double deltp = dv[ibn].press - param.PsOutlet;
      double rhob = cv[ibn].dens - deltp / (dv[ibn].cs * dv[ibn].cs);
      double ub = u + sxn * deltp * rrhoc;
      double vb = v + syn * deltp * rrhoc;
      double vnd = ub * sxn + vb * syn;

      if (vnd < 0.0) {
        ub = (((u) < (0)) ? (-1.0) : (1.0)) *
             std::max(std::abs(ub), std::abs(u));
        vb = (((v) < (0)) ? (-1.0) : (1.0)) *
             std::max(std::abs(vb), std::abs(v));
      }
      cv[idn].dens = rhob;
      cv[idn].xmom = rhob * ub;
      cv[idn].ymom = rhob * vb;
      cv[idn].ener = param.PsOutlet / gam1 + 0.5 * rhob * (ub * ub + vb * vb);
    } else {
      cv[idn].dens = cv[ibn].dens;
      cv[idn].xmom = cv[ibn].xmom;
      cv[idn].ymom = cv[ibn].ymom;
      cv[idn].ener = cv[ibn].ener;
    }
    ConvToDepend(idn);
  }
}

void FVMSolver::BoundWallVisc(int beg, int end) {
  for (int ib = beg; ib <= end; ++ib) {
    int ibn = geom.boundNode[ib].node;

    cv[ibn].xmom = 0.0;
    cv[ibn].ymom = 0.0;

    ConvToDepend(ibn);
  }
}

void FVMSolver::Gradients() {
  double fcx[4];
  double fcy[4];

  PRIM_VAR init;
  init.dens = 0.0;
  init.velx = 0.0;
  init.vely = 0.0;
  init.press = 0.0;

  gradx.assign(geom.totNodes, init);
  grady.assign(geom.totNodes, init);

  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;

    double rav = 0.5 * (cv[i].dens + cv[j].dens);
    double uav = 0.5 * (cv[i].xmom / cv[i].dens + cv[j].xmom / cv[j].dens);
    double vav = 0.5 * (cv[i].ymom / cv[i].dens + cv[j].ymom / cv[j].dens);
    double pav = 0.5 * (dv[i].press + dv[j].press);

    double sx = geom.sij[ie].x;
    double sy = geom.sij[ie].y;

    fcx[0] = rav * sx;
    fcx[1] = uav * sx;
    fcx[2] = vav * sx;
    fcx[3] = pav * sx;

    fcy[0] = rav * sy;
    fcy[1] = uav * sy;
    fcy[2] = vav * sy;
    fcy[3] = pav * sy;

    gradx[i].dens += fcx[0];
    gradx[i].velx += fcx[1];
    gradx[i].vely += fcx[2];
    gradx[i].press += fcx[3];
    gradx[j].dens -= fcx[0];
    gradx[j].velx -= fcx[1];
    gradx[j].vely -= fcx[2];
    gradx[j].press -= fcx[3];

    grady[i].dens += fcy[0];
    grady[i].velx += fcy[1];
    grady[i].vely += fcy[2];
    grady[i].press += fcy[3];
    grady[j].dens -= fcy[0];
    grady[j].velx -= fcy[1];
    grady[j].vely -= fcy[2];
    grady[j].press -= fcy[3];
  }

  int ibegf = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendf = geom.ibound[ib].bfaceIndex;
    bool flag = true;
    if (geom.BoundTypes[ib] >= 500 && geom.BoundTypes[ib] < 600) {
      flag = false;
    }
    if (geom.BoundTypes[ib] >= 700 && geom.BoundTypes[ib] < 800) {
      flag = false;
    }

    if (flag) {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        int i = geom.boundFace[ibf].nodei;
        int j = geom.boundFace[ibf].nodej;
        double sx = geom.sbf[ibf].x / 12.0;
        double sy = geom.sbf[ibf].y / 12.0;

        double rav = 5.0 * cv[i].dens + cv[j].dens;
        double uav = 5.0 * cv[i].xmom / cv[i].dens + cv[j].xmom / cv[j].dens;
        double vav = 5.0 * cv[i].ymom / cv[i].dens + cv[j].ymom / cv[j].dens;
        double pav = 5.0 * dv[i].press + dv[j].press;

        gradx[i].dens += rav * sx;
        gradx[i].velx += uav * sx;
        gradx[i].vely += vav * sx;
        gradx[i].press += pav * sx;

        grady[i].dens += rav * sy;
        grady[i].velx += uav * sy;
        grady[i].vely += vav * sy;
        grady[i].press += pav * sy;

        rav = 5.0 * cv[j].dens + cv[i].dens;
        uav = 5.0 * cv[j].xmom / cv[j].dens + cv[i].xmom / cv[i].dens;
        vav = 5.0 * cv[j].ymom / cv[j].dens + cv[i].ymom / cv[i].dens;
        pav = 5.0 * dv[j].press + dv[i].press;

        gradx[j].dens += rav * sx;
        gradx[j].velx += uav * sx;
        gradx[j].vely += vav * sx;
        gradx[j].press += pav * sx;

        grady[j].dens += rav * sy;
        grady[j].velx += uav * sy;
        grady[j].vely += vav * sy;
        grady[j].press += pav * sy;
      }
    }
    ibegf = iendf + 1;
  }

  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    if (geom.BoundTypes[ib] >= 500 && geom.BoundTypes[ib] < 600) {
      if ((geom.BoundTypes[ib] - 500) < 2) {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.boundNode[ibn].node;
          gradx[i].dens = 0.0;
          grady[i].velx = 0.0;
          gradx[i].vely = 0.0;
          gradx[i].press = 0.0;
        }
      } else {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.boundNode[ibn].node;
          grady[i].dens = 0.0;
          grady[i].velx = 0.0;
          gradx[i].vely = 0.0;
          grady[i].press = 0.0;
        }
      }
    }
    ibegn = iendn + 1;
  }

  PeriodicPrim(gradx);
  PeriodicPrim(grady);

  for (int i = 0; i < geom.phyNodes; ++i) {
    gradx[i].dens /= geom.vol[i];
    gradx[i].velx /= geom.vol[i];
    gradx[i].vely /= geom.vol[i];
    gradx[i].press /= geom.vol[i];

    grady[i].dens /= geom.vol[i];
    grady[i].velx /= geom.vol[i];
    grady[i].vely /= geom.vol[i];
    grady[i].press /= geom.vol[i];
  }
}

void FVMSolver::PeriodicPrim(std::vector<PRIM_VAR> &var) {
  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    if (geom.BoundTypes[ib] >= 700 && geom.BoundTypes[ib] < 800) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.boundNode[ibn].node;
        int j = geom.boundNode[ibn].dummy;

        var[i].dens += var[j].dens;
        var[i].velx += var[j].velx;
        var[i].vely += var[j].vely;
        var[i].press += var[j].press;
        var[j].dens = var[i].dens;
        var[j].velx = var[i].velx;
        var[j].vely = var[i].vely;
        var[j].press = var[i].press;
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::PeriodicCons(std::vector<CONS_VAR> &var) {
  int ibegn = 0.0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    if (geom.BoundTypes[ib] >= 700 && geom.BoundTypes[ib] < 800) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.boundNode[ibn].node;
        int j = geom.boundNode[ibn].dummy;
        var[i].dens += var[j].dens;
        var[i].xmom += var[j].xmom;
        var[i].ymom += var[j].ymom;
        var[i].ener += var[j].ener;
        var[j].dens = var[i].dens;
        var[j].xmom = var[i].xmom;
        var[j].ymom = var[i].ymom;
        var[j].ener = var[i].ener;
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::PeriodicInt(std::vector<int> &var) {
  int i, j, ib, ibn, ibegn, iendn;

  ibegn = 0;

  for (ib = 0; ib < geom.numBoundSegs; ib++) {
    iendn = geom.ibound[ib].bnodeIndex;
    if (geom.BoundTypes[ib] >= 700 && geom.BoundTypes[ib] < 800) {
      for (ibn = ibegn; ibn <= iendn; ibn++) {
        i = geom.boundNode[ibn].node;
        j = geom.boundNode[ibn].dummy;
        var[i] += var[j];
        var[j] = var[i];
      }
    }
    ibegn = iendn + 1;
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

void FVMSolver::Timestep() {
  if (param.equationtype_ == preprocess::equationType::Euler) {
    for (int i = 0; i < geom.phyNodes; ++i) {
      double sx = geom.sproj[i].x;
      double sy = geom.sproj[i].y;
      double u = std::abs(cv[i].xmom / cv[i].dens);
      double v = std::abs(cv[i].ymom / cv[i].dens);
      double vc = sx * u + sy * v;
      double cs = dv[i].cs * (sx + sy);
      timeSteps[i] = geom.vol[i] / (vc + cs);
    }
  } else {
    double cfac = 2.0;
    for (int i = 0; i < geom.phyNodes; ++i) {
      double sx = geom.sproj[i].x;
      double sy = geom.sproj[i].y;
      double ds = sx + sy;
      double rrho = 1.0 / cv[i].dens;
      double u = std::abs(cv[i].xmom / cv[i].dens);
      double v = std::abs(cv[i].ymom / cv[i].dens);
      double fmue = dvlam[i].mu / param.Prandtl;
      double f1 = (4.0 * rrho) / 3.0;
      double f2 = dv[i].gamma * rrho;
      double fac = std::max(f1, f2);
      double dtv = (fac * fmue) / geom.vol[i];
      double vc = sx * u + sy * v;
      double cs = dv[i].cs * (sx + sy);
      double lambdac = vc + cs;
      double lambdav = dtv * (sx * sx + sy * sy);
      timeSteps[i] = geom.vol[i] / (lambdac + cfac * lambdav);
    }
  }

  if (param.timestep_ == preprocess::timeStep::global) {
    double tsmin = 1e+33;
    for (int i = 0; i < geom.phyNodes; ++i) {
      tsmin = std::min(tsmin, timeSteps[i]);
    }
    for (int i = 0; i < geom.phyNodes; ++i) {
      timeSteps[i] = tsmin;
    }
  }
}

void FVMSolver::solve() {
  for (int i = 0; i < geom.totNodes; ++i) {
    cvOld[i] = cv[i];
  }

  Timestep();

  for (int irk = 0; irk < param.temperalStages; ++irk) {
    if (param.dissipationEval[irk]) {
      DissipInit(irk, param.dissipationBlend[irk]);
    }

    if (param.dissipationEval[irk]) {
      if (param.equationtype_ == preprocess::equationType::Euler) {
        Gradients();
      }
      limiter->limiterInit(param, geom, cv, dv, umin, umax, lim);
      limiter->limiterUpdate(param, geom, cv, dv, umin, umax, lim, gradx,
                             grady);

      DissipRoe2(param.dissipationBlend[irk]);
    }

    FluxRoe2();
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

    if (param.imResiSmooth > 0.0) {
      Irsmoo();
      ZeroRes();
    }

    for (int i = 0; i < geom.phyNodes; ++i) {
      cv[i].dens = cvOld[i].dens - rhs[i].dens;
      cv[i].xmom = cvOld[i].xmom - rhs[i].xmom;
      cv[i].ymom = cvOld[i].ymom - rhs[i].ymom;
      cv[i].ener = cvOld[i].ener - rhs[i].ener;
    }

    ConvToDependAll();

    BoundaryConditions();
  }
}

void FVMSolver::ZeroRes() {
  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    if (geom.BoundTypes[ib] >= 500 && geom.BoundTypes[ib] < 600) {
      if (geom.BoundTypes[ib] - 500 < 2) {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.boundNode[ibn].node;
          rhs[i].xmom = 0.0;
        }
      } else {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.boundNode[ibn].node;
          rhs[i].ymom = 0.0;
        }
      }
    } else if ((geom.BoundTypes[ib] >= 300 && geom.BoundTypes[ib] < 400) &&
               param.equationtype_ == preprocess::equationType::NavierStokes) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.boundNode[ib].node;
        rhs[i].xmom = 0.0;
        rhs[i].ymom = 0.0;
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::Irsmoo() {
  for (int i = 0; i < geom.totNodes; ++i) {
    nContr[i] = 0;
    rhsOld[i] = rhs[i];
  }

  for (int itirs = 1; itirs <= param.numOfIterSmooth; ++itirs) {
    for (int i = 0; i < geom.phyNodes; ++i) {
      rhsIter[i].dens = 0.0;
      rhsIter[i].xmom = 0.0;
      rhsIter[i].ymom = 0.0;
      rhsIter[i].ener = 0.0;
    }
    if (itirs == 1) {
      for (int ie = 0; ie < geom.phyEdges; ++ie) {
        int i = geom.edge[ie].nodei;
        int j = geom.edge[ie].nodej;
        nContr[i]++;
        nContr[j]++;

        rhsIter[i].dens += rhs[j].dens;
        rhsIter[i].xmom += rhs[j].xmom;
        rhsIter[i].ymom += rhs[j].ymom;
        rhsIter[i].ener += rhs[j].ener;

        rhsIter[j].dens += rhs[i].dens;
        rhsIter[j].xmom += rhs[i].xmom;
        rhsIter[j].ymom += rhs[i].ymom;
        rhsIter[j].ener += rhs[i].ener;
      }

      PeriodicInt(nContr);
    } else {
      for (int ie = 0; ie < geom.phyEdges; ++ie) {
        int i = geom.edge[ie].nodei;
        int j = geom.edge[ie].nodej;
        rhsIter[i].dens += rhs[j].dens;
        rhsIter[i].xmom += rhs[j].xmom;
        rhsIter[i].ymom += rhs[j].ymom;
        rhsIter[i].ener += rhs[j].ener;

        rhsIter[j].dens += rhs[i].dens;
        rhsIter[j].xmom += rhs[i].xmom;
        rhsIter[j].ymom += rhs[i].ymom;
        rhsIter[j].ener += rhs[i].ener;
      }
    }
    PeriodicCons(rhsIter);

    for (int i = 0; i < geom.phyNodes; ++i) {
      double den = 1.0 / (1.0 + ((double)nContr[i]) * param.imResiSmooth);
      rhs[i].dens =
          (rhsIter[i].dens * param.imResiSmooth + rhsOld[i].dens) * den;
      rhs[i].xmom =
          (rhsIter[i].xmom * param.imResiSmooth + rhsOld[i].xmom) * den;
      rhs[i].ymom =
          (rhsIter[i].ymom * param.imResiSmooth + rhsOld[i].ymom) * den;
      rhs[i].ener =
          (rhsIter[i].ener * param.imResiSmooth + rhsOld[i].ener) * den;
    }
  }
}

void FVMSolver::DensityChange(double &drho, double &drmax, int &idrmax) {
  int i;
  double dr;

  drho = 0.0;
  drmax = 0.0;
  idrmax = 0;
  for (i = 0; i < geom.phyNodes; i++) {
    dr = cv[i].dens - cvOld[i].dens;
    drho += dr * dr;
    if (std::abs(dr) >= drmax) {
      drmax = std::abs(dr);
      idrmax = i;
    }
  }

  drho = std::sqrt(drho);
}

void FVMSolver::Convergence() {
  int idr;
  double drmax;

  DensityChange(drho, drmax, idr);

  if (iter == 1) {
    drho1 = drho;
    drho1 = std::max(drho1, 1.0e-33);
    drho = 1.0;
  } else {
    drho = drho / drho1;
    drho = std::max(drho, 1.0e-33);
  }

  if (param.flowtype_ == preprocess::flowType::External) {
    Forces();
  } else {
    MassFlow();
  }
  if (param.flowtype_ == preprocess::flowType::External) {
    std::cout << "iter: " << iter << " " << "drho: " << std::log10(drho) << " "
              << "drmax: " << drmax << " " << "idr: " << idr << " "
              << "cl: " << cl << " " << "cd: " << cd << " " << "cm: " << cm
              << "\n";
  } else {
    std::cout << "iter: " << iter << " " << "drho: " << std::log10(drho) << " "
              << "drmax: " << drmax << " " << "idr: " << idr << " "
              << "Mass Flow: " << mflow << " " << "Mass Flow Ratio: " << mfratio
              << "\n";
  }
}

void FVMSolver::writeTecplotDat() {
  std::ofstream outFile(param.title + ".dat");
  outFile << "TITLE = \"Unstructured Triangle Mesh\"\n";
  outFile
      << "VARIABLES = \"X\", \"Y\", \"Density\", \"U\", \"V\", \"P\", \"T\"\n";
  outFile << "ZONE N=" << geom.phyNodes << ", E=" << geom.numTria
          << ", F=FEPOINT, ET=TRIANGLE\n";

  for (int i = 0; i < geom.phyNodes; ++i) {
    outFile << geom.coords[i].x << " " << geom.coords[i].y << " " << cv[i].dens
            << " " << cv[i].xmom / cv[i].dens << " " << cv[i].ymom / cv[i].dens
            << " " << dv[i].press << " " << dv[i].temp << "\n";
  }

  for (int i = 0; i < geom.numTria; ++i) {
    outFile << geom.tria[i].node[0] + 1 << " " << geom.tria[i].node[1] + 1
            << " " << geom.tria[i].node[2] + 1 << "\n";
  }
  outFile << "\n";
  outFile.close();
}

void FVMSolver::Forces() {
  double cx = 0.0;
  double cy = 0.0;
  cm = 0.0;
  int ibegf = 0.0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendf = geom.ibound[ib].bfaceIndex;
    if (geom.BoundTypes[ib] >= 300 && geom.BoundTypes[ib] < 500) {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        int n1 = geom.boundFace[ibf].nodei;
        int n2 = geom.boundFace[ibf].nodej;
        double sx = geom.sbf[ibf].x;
        double sy = geom.sbf[ibf].y;
        double pwall = 0.5 * (dv[n1].press + dv[n2].press);
        double cp = 2.0 * (pwall - param.PsInfinity) /
                    (param.RhoInfinity * param.velInfinity * param.velInfinity);
        double xa =
            (0.5 * (geom.coords[n1].x + geom.coords[n2].x) - param.xRefPoint) /
            param.chord;
        double ya =
            (0.5 * (geom.coords[n1].y + geom.coords[n2].y) - param.yRefPoint) /
            param.chord;
        double dcy = sy * cp;
        double dcx = sx * cp;
        cy = cy + dcy;
        cx = cx + dcx;
        cm = cm + dcx * ya - dcy * xa;
      }
    }
    ibegf = iendf + 1;
  }

  cl = cy * std::cos(param.Aoa) - cx * std::sin(param.Aoa);
  cd = cy * std::sin(param.Aoa) + cx * std::cos(param.Aoa);
}

void FVMSolver::MassFlow() {
  double massin = 0.0;
  double massout = 0.0;
  double mass = 0.0;
  mflow = 0.0;
  mfratio = 0.0;
  bool in = false;
  bool out = false;

  int ibegf = 0.0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendf = geom.ibound[ib].bfaceIndex;
    if (geom.BoundTypes[ib] >= 100 && geom.BoundTypes[ib] < 300) {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        int n1 = geom.boundFace[ibf].nodei;
        int n2 = geom.boundFace[ibf].nodej;
        double sx = geom.sbf[ibf].x;
        double sy = geom.sbf[ibf].y;
        mass = 0.5 * ((cv[n1].xmom + cv[n2].xmom) * sx +
                      (cv[n1].ymom + cv[n2].ymom) * sy);
        if (geom.BoundTypes[ib] >= 100 && geom.BoundTypes[ib] < 200) {
          massin = massin - mass;
          in = true;
        } else {
          massout = massout + mass;
          out = true;
        }
      }
    }
    ibegf = iendf + 1;
  }

  if (in) {
    mflow = massin;
    mfratio = massout / massin;
  }
  if (!in && out) {
    mflow = massout;
    mfratio = 1.0;
  }
}
}  // namespace solver