#include <solver/FVMSolver.h>

#include <cmath>
#include <iostream>

#include "pre/parameter.h"

namespace solver {
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

void FVMSolver::Forces() {
  double cx = 0.0;
  double cy = 0.0;
  cm = 0.0;
  int ibegf = 0.0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendf = geom.ibound[ib].bfaceIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if ((type == preprocess::BoundaryType::EulerWall) ||
        (type == preprocess::BoundaryType::NoSlipWall)) {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        int n1 = geom.boundaryFace[ibf].nodei;
        int n2 = geom.boundaryFace[ibf].nodej;
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
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Inflow ||
        type == preprocess::BoundaryType::Outflow) {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        int n1 = geom.boundaryFace[ibf].nodei;
        int n2 = geom.boundaryFace[ibf].nodej;
        double sx = geom.sbf[ibf].x;
        double sy = geom.sbf[ibf].y;
        mass = 0.5 * ((cv[n1].xmom + cv[n2].xmom) * sx +
                      (cv[n1].ymom + cv[n2].ymom) * sy);
        if (type == preprocess::BoundaryType::Inflow) {
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