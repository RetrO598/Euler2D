#include <solver/FVMSolver.h>

#include <cmath>

namespace solver {
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
}  // namespace solver