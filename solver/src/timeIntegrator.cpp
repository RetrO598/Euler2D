#include "solver/FVMSolver.h"
#include "solver/timeIntegrator.hpp"
namespace solver {

RungeKuttaTimeIntegrator::RungeKuttaTimeIntegrator(FVMSolver& solver)
    : solver(solver) {}

void RungeKuttaTimeIntegrator::timeAdvance() {
  solver.assignCVold();
  solver.Timestep();

  for (int irk = 0; irk < solver.param.temperalStages; ++irk) {
    solver.updateResidualRK(irk);

    if (solver.param.imResiSmooth > 0.0) {
      solver.Irsmoo();
      solver.ZeroRes();
    }

    solver.updateCV();
    solver.ConvToDependAll();

    solver.WallVisc();
  }
}

LUSGSIntegrator::LUSGSIntegrator(FVMSolver& solver) : solver(solver) {}

void LUSGSIntegrator::timeAdvance() {
  solver.assignCVold();
  solver.Timestep();
  solver.computeWaveSpeed();
  solver.computeJacobianDiag(solver.param.lusgsParameter);
  solver.computeResidualLUSGS();
  if (solver.param.imResiSmooth > 0.0) {
    solver.Irsmoo();
    solver.ZeroRes();
  }
  solver.lowerSweep();
  solver.upperSweep();
  solver.LUSGSupdate();
  solver.ConvToDependAll();
  solver.WallVisc();
}
}  // namespace solver