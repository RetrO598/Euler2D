#pragma once

#include <pre/geometry.h>
#include <pre/parameter.h>
#include <solver/variableDef.h>

#include <vector>
namespace solver {
class BaseNumeric {
public:
  BaseNumeric(const preprocess::parameter &param,
              const preprocess::Geometry &geom, std::vector<CONS_VAR> &cv,
              std::vector<DEPEND_VAR> &dv, std::vector<CONS_VAR> &diss,
              std::vector<CONS_VAR> &rhs, std::vector<PRIM_VAR> &lim,
              std::vector<PRIM_VAR> &gradx, std::vector<PRIM_VAR> &grady);

  void DissipInit(const int &irk, const double &beta);

  void FluxWalls();

  virtual void DissipNumeric(const double &beta){};

  virtual void FluxNumeric(){};

protected:
  const preprocess::parameter &param;
  const preprocess::Geometry &geom;
  std::vector<CONS_VAR> &cv;
  std::vector<DEPEND_VAR> &dv;
  std::vector<CONS_VAR> &diss;
  std::vector<CONS_VAR> &rhs;
  std::vector<PRIM_VAR> &lim;
  std::vector<PRIM_VAR> &gradx;
  std::vector<PRIM_VAR> &grady;

  friend class FVMSolver;
};

class NumericRoe : public BaseNumeric {
public:
  NumericRoe(const preprocess::parameter &param,
             const preprocess::Geometry &geom, std::vector<CONS_VAR> &cv,
             std::vector<DEPEND_VAR> &dv, std::vector<CONS_VAR> &diss,
             std::vector<CONS_VAR> &rhs, std::vector<PRIM_VAR> &lim,
             std::vector<PRIM_VAR> &gradx, std::vector<PRIM_VAR> &grady);

  void DissipNumeric(const double &beta) override;

  void FluxNumeric() override;

  double EntropyCorr(const double &z, const double &d);
};

class NumericSLAU2 : public BaseNumeric {
public:
  NumericSLAU2(const preprocess::parameter &param,
               const preprocess::Geometry &geom, std::vector<CONS_VAR> &cv,
               std::vector<DEPEND_VAR> &dv, std::vector<CONS_VAR> &diss,
               std::vector<CONS_VAR> &rhs, std::vector<PRIM_VAR> &lim,
               std::vector<PRIM_VAR> &gradx, std::vector<PRIM_VAR> &grady);

  void DissipNumeric(const double &beta) override;
  void FluxNumeric() override;

  inline double pressureFuncLeft(const double &m) {
    if (std::abs(m) >= 1.0) {
      return 0.5 * (1 + (m >= 0 ? 1 : -1));
    } else {
      return 0.25 * (m + 1) * (m + 1) * (2 - m);
    }
  }
  inline double pressureFuncRight(const double &m) {
    if (std::abs(m) >= 1.0) {
      return 0.5 * (1 - (m >= 0 ? 1 : -1));
    } else {
      return 0.25 * (m - 1) * (m - 1) * (2 + m);
    }
  }
};

class NumericAUSM : public BaseNumeric {
public:
  NumericAUSM(const preprocess::parameter &param,
              const preprocess::Geometry &geom, std::vector<CONS_VAR> &cv,
              std::vector<DEPEND_VAR> &dv, std::vector<CONS_VAR> &diss,
              std::vector<CONS_VAR> &rhs, std::vector<PRIM_VAR> &lim,
              std::vector<PRIM_VAR> &gradx, std::vector<PRIM_VAR> &grady);
  void DissipNumeric(const double &beta) override;
  void FluxNumeric() override;
};
} // namespace solver