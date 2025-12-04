#pragma once

#include "solver/FVMSolver.h"
namespace solver {
class TimeIntegrator {
 public:
  virtual void timeAdvance() = 0;
  virtual ~TimeIntegrator() = default;
};

class RungeKuttaTimeIntegrator : public TimeIntegrator {
 public:
  RungeKuttaTimeIntegrator(FVMSolver& solver);
  void timeAdvance() override;
  ~RungeKuttaTimeIntegrator() = default;

 private:
  FVMSolver& solver;
};

}  // namespace solver