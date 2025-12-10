#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "pre/parameter.h"
#include "solver/FVMSolver.h"
#include "solver/turbTemplate.hpp"
namespace solver {
void FVMSolver::TurbGradients() {
  gradTurbX.assign(geom.phyNodes, 0.0);
  gradTurbY.assign(geom.phyNodes, 0.0);

  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;

    double sx = geom.sij[ie].x;
    double sy = geom.sij[ie].y;

    double ave = 0.5 * (turbVar[i].nu_turb + turbVar[j].nu_turb);
    gradTurbX[i] += ave * sx;
    gradTurbX[j] -= ave * sx;

    gradTurbY[i] += ave * sy;
    gradTurbY[j] -= ave * sy;

    if (std::isnan(ave)) {
      std::cout << "not a number from gradient\n";
      exit(1);
    }
  }

  int ibegf = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendf = geom.ibound[ib].bfaceIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    bool flag = true;
    if (type == preprocess::BoundaryType::Symmetric) {
      flag = false;
    }
    if (type == preprocess::BoundaryType::Periodic) {
      flag = false;
    }

    if (flag) {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        int i = geom.boundaryFace[ibf].nodei;
        int j = geom.boundaryFace[ibf].nodej;
        double sx = geom.sbf[ibf].x / 12.0;
        double sy = geom.sbf[ibf].y / 12.0;

        double ave = 5.0 * turbVar[i].nu_turb + turbVar[j].nu_turb;

        gradTurbX[i] += ave * sx;
        gradTurbY[i] += ave * sy;

        ave = 5.0 * turbVar[j].nu_turb + turbVar[i].nu_turb;

        gradTurbX[j] += ave * sx;
        gradTurbY[j] += ave * sy;

        if (std::isnan(ave)) {
          std::cout << "not a number from gradient\n";
          exit(1);
        }
      }
    }
    ibegf = iendf + 1;
  }

  ibegf = 0;
  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    int iendf = geom.ibound[ib].bfaceIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Symmetric) {
      double sx = 0.0;
      double sy = 0.0;
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        sx += geom.sbf[ibf].x;
        sy += geom.sbf[ibf].y;
      }
      if (std::abs(sx) > std::abs(sy)) {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.vertexList[ibn].nodeIdx;
          gradTurbX[i] = 0.0;
        }
      } else {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.vertexList[ibn].nodeIdx;
          gradTurbY[i] = 0.0;
        }
      }
    }
    ibegn = iendn + 1;
    ibegf = iendf + 1;
  }

  PeriodicVisc(gradTurbX, gradTurbY);

  for (int i = 0; i < geom.phyNodes; ++i) {
    gradTurbX[i] /= geom.vol[i];
    gradTurbY[i] /= geom.vol[i];
  }
}

void FVMSolver::TurbConvection() {
  for (int i = 0; i < geom.phyNodes; ++i) {
    rhsTurb[i].nu_turb = -dissTurb[i].nu_turb;
  }

  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;

    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);

    double nx = geom.sij[ie].x / ds;
    double ny = geom.sij[ie].y / ds;

    UpwindTurbSA::ComputeResidual(nx, ny, ds, turbVar[i], turbVar[j], cv[i],
                                  cv[j], rhsTurb[i], rhsTurb[j]);
  }
}

void FVMSolver::TurbViscous() {
  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;

    double sx = geom.sij[ie].x;
    double sy = geom.sij[ie].y;

    double nuLami = dvlam[i].mu / cv[i].dens;
    double nuLamj = dvlam[j].mu / cv[j].dens;

    double nuLamAve = 0.5 * (nuLamj + nuLami);

    double nuTurbAve = 0.5 * (turbVar[i].nu_turb + turbVar[j].nu_turb);

    double gradxAve = 0.5 * (gradTurbX[i] + gradTurbX[j]);
    double gradyAve = 0.5 * (gradTurbY[i] + gradTurbY[j]);

    double sigma = 2.0 / 3.0;

    double fv;
    fv = sx * (nuTurbAve + nuLamAve) * gradxAve / sigma +
         sy * (nuTurbAve + nuLamAve) * gradyAve / sigma;

    dissTurb[i].nu_turb += fv;
    dissTurb[j].nu_turb -= fv;

    if (std::isnan(dvlam[i].mu) || std::isnan(dvlam[j].mu)) {
      std::cout << "not a number from viscous" << "\n";
      std::cout << i << " " << dvlam[i].mu << "\n";
      std::cout << j << " " << dvlam[j].mu << "\n";
      exit(1);
    }
  }
}

void FVMSolver::TurbSource() {
  for (int i = 0; i < geom.phyNodes; ++i) {
    double cv1_3 = 357.911;
    double k2 = 0.1681;
    double cb1 = 0.1355;
    double cw2 = 0.3;
    double cw3_6 = 64.;
    double sigma = 2. / 3.;
    double cb2 = 0.622;
    double cb2_sigma = cb2 / sigma;
    double cw1 = cb1 / k2 + (1 + cb2) / sigma;

    double density = cv[i].dens;
    double nuLam = dvlam[i].mu / density;

    double production = 0.0;
    double destruction = 0.0;
    double crossProduction = 0.0;

    double vorticity =
        (gradx[i].vely - grady[i].velx) * (gradx[i].vely - grady[i].velx);
    double omega = std::sqrt(vorticity);

    if (geom.pointList[i].getWallDistance() > 1e-10) {
      double dist = geom.pointList[i].getWallDistance();
      double dist2 = dist * dist;
      double Ji = turbVar[i].nu_turb / nuLam;
      double Ji2 = Ji * Ji;
      double Ji3 = Ji2 * Ji;
      double fv1 = Ji3 / (Ji3 + cv1_3);
      double fv2 = 1.0 - Ji / (1.0 + Ji * fv1);
      double S = omega;
      double inv_k2_d2 = 1.0 / (k2 * dist2);
      double Shat = S + turbVar[i].nu_turb * fv2 * inv_k2_d2;
      double inv_Shat = 1.0 / std::max(Shat, 1.0e-10);
      production = cb1 * Shat * turbVar[i].nu_turb;

      double r = std::min(turbVar[i].nu_turb * inv_Shat * inv_k2_d2, 10.0);
      double g = r + cw2 * (r * r * r * r * r * r - r);
      double g6 = g * g * g * g * g * g;
      double glim = std::pow((1.0 + cw3_6) / (g6 + cw3_6), 1.0 / 6.0);
      double fw = g * glim;

      destruction = cw1 * fw * turbVar[i].nu_turb * turbVar[i].nu_turb / dist2;

      double norm2_Grad = 0.0;
      norm2_Grad = (gradTurbX[i] * gradTurbX[i] + gradTurbY[i] * gradTurbY[i]);

      crossProduction = cb2_sigma * norm2_Grad;

      rhsTurb[i].nu_turb -=
          (production - destruction + crossProduction) * geom.vol[i];

      if (std::isnan(production) || std::isnan(destruction) ||
          std::isnan(crossProduction)) {
        std::cout << "not a number from source" << "\n";
        exit(1);
      }
    }
  }
}

void FVMSolver::TurbDissInit() {
  for (int i = 0; i < geom.phyNodes; ++i) {
    dissTurb[i].nu_turb = 0.0;
  }
}

void FVMSolver::initTurbSolver() {
  if (param.equationtype_ == preprocess::equationType::RANS) {
    for (int i = 0; i < geom.phyNodes; ++i) {
      turbVar[i].nu_turb = 0.1 * dvlam[i].mu / cv[i].dens;
    }
  }
}
}  // namespace solver