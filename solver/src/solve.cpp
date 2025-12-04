#include <solver/FVMSolver.h>

#include <cmath>

#include "pre/parameter.h"

namespace solver {
// void FVMSolver::solve() {
//   for (int i = 0; i < geom.phyNodes; ++i) {
//     cvOld[i] = cv[i];
//   }

//   Timestep();

//   for (int irk = 0; irk < param.temperalStages; ++irk) {
//     numeric->DissipInit();

//     if (param.equationtype_ == preprocess::equationType::NavierStokes) {
//       GradientsVisc();
//       fluxVisc();
//     } else if (param.equationtype_ == preprocess::equationType::Euler) {
//       Gradients();
//     }

//     limiter->limiterInit();
//     limiter->limiterUpdate();

//     numeric->FluxNumeric();
//     BoundaryConditions();
//     ZeroRes();
//     PeriodicCons(rhs);

//     double fac = param.stageCoeff[irk] * param.CFL;
//     for (int i = 0; i < geom.phyNodes; ++i) {
//       double adtv = fac * timeSteps[i] / geom.vol[i];
//       rhs[i].dens *= adtv;
//       rhs[i].xmom *= adtv;
//       rhs[i].ymom *= adtv;
//       rhs[i].ener *= adtv;
//     }

//     if (param.imResiSmooth > 0.0) {
//       Irsmoo();
//       ZeroRes();
//     }

//     for (int i = 0; i < geom.phyNodes; ++i) {
//       cv[i].dens = cvOld[i].dens - rhs[i].dens;
//       cv[i].xmom = cvOld[i].xmom - rhs[i].xmom;
//       cv[i].ymom = cvOld[i].ymom - rhs[i].ymom;
//       cv[i].ener = cvOld[i].ener - rhs[i].ener;
//     }

//     ConvToDependAll();

//     WallVisc();
//   }
// }
}  // namespace solver