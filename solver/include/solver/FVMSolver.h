#pragma once

#include <pre/geometry.h>
#include <pre/parameter.h>
#include <solver/limiter.h>
#include <solver/numeric.h>
#include <solver/variableDef.h>

#include <vector>
namespace solver {

class FVMSolver {
public:
  FVMSolver(preprocess::parameter &parameter,
            const preprocess::Geometry &geometry);

  void initSolver();
  void readSolution();
  void ConvToDependAll();
  void ConvToDepend(int i);

  void BoundaryConditions();
  void BoundInflow(int beg, int end);
  void BoundOutflow(int beg, int end);
  void BoundFarfield(int beg, int end);
  void BoundWallVisc(int beg, int end);

  void Gradients();

  void PeriodicCons(std::vector<CONS_VAR> &var);
  void PeriodicPrim(std::vector<PRIM_VAR> &var);
  void PeriodicInt(std::vector<int> &var);
  void PeriodicVisc();

  void solve();

  void Timestep();

  void ZeroRes();

  void Irsmoo();

  void DensityChange(double &drho, double &drmax, int &idrmax);

  void Convergence();

  void Forces();
  void MassFlow();

  inline bool Converged() {
    if (iter >= param.maxIteration || drho <= param.convTol) {
      return true;
    }
    return false;
  }

  void writeTecplotDat();
  void writeLineDat();

  preprocess::parameter &param;
  const preprocess::Geometry &geom;
  BaseLimiter *limiter;
  BaseNumeric *numeric;

  int iter;

private:
  std::vector<CONS_VAR> cv;
  std::vector<CONS_VAR> cvOld;
  std::vector<CONS_VAR> diss;
  std::vector<CONS_VAR> rhs;

  std::vector<PRIM_VAR> lim;
  std::vector<PRIM_VAR> gradx;
  std::vector<PRIM_VAR> grady;
  std::vector<PRIM_VAR> umin;
  std::vector<PRIM_VAR> umax;

  std::vector<double> gradTx;
  std::vector<double> gradTy;

  std::vector<DEPEND_VAR> dv;

  std::vector<VISC_VAR> dvlam;

  std::vector<double> timeSteps;

  std::vector<CONS_VAR> rhsIter;
  std::vector<CONS_VAR> rhsOld;
  std::vector<int> nContr;

  double cl;
  double cd;
  double cm;
  double mflow;
  double mfratio;
  double drho;
  double drho1;
};
} // namespace solver