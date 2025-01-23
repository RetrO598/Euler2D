#include <solver/FVMSolver.h>

#include <cmath>

namespace solver {
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
      limiter->limiterInit();
      limiter->limiterUpdate();

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
}  // namespace solver