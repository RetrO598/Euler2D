#pragma once

#include <limits>

#include "solver/FVMSolver.h"
#include "solver/variableDef.h"
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

class LUSGSIntegrator : public TimeIntegrator {
 public:
  LUSGSIntegrator(FVMSolver& solver);
  void timeAdvance() override;
  ~LUSGSIntegrator() = default;

 private:
  FVMSolver& solver;
};

struct fluxDifference {
  inline static CONS_VAR computeDifference(CONS_VAR& u, double nx, double ny,
                                           double gamma) {
    CONS_VAR result;
    double rho = u.dens;
    double et = u.ener;
    double temp = 0.0;
    double contrav = 0.0;

    contrav = (u.xmom / rho * nx + u.ymom / rho * ny);
    temp = (u.xmom * u.xmom + u.ymom * u.ymom);

    double p = (gamma - 1.0) * (et - 0.5 * temp / rho);
    if (p < std::numeric_limits<double>::epsilon()) {
      p = std::numeric_limits<double>::epsilon();
      et = p / (gamma - 1.0) + 0.5 * temp / rho;
      u.ener = et;
    }

    double ht = et + p;
    result.dens = rho * contrav;
    result.xmom = u.xmom * contrav + p * nx;
    result.ymom = u.ymom * contrav + p * ny;
    result.ener = ht * contrav;

    return result;
  }
};
}  // namespace solver