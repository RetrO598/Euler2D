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

  interfaceVar Interpolate(const int &i, const int &j);

  virtual void DissipNumeric(const double &beta) {};

  virtual void FluxNumeric() {};

  virtual void ComputeResidual(const preprocess::parameter &param, double nx,
                               double ny, double ds, PRIM_VAR vari,
                               PRIM_VAR varj, double gammai, double gammaj,
                               CONS_VAR &resi, CONS_VAR &resj) {};

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
  using BaseNumeric::BaseNumeric;

  void DissipNumeric(const double &beta) override;

  void FluxNumeric() override;

  double EntropyCorr(const double &z, const double &d);
  void ComputeResidual(const preprocess::parameter &param, double nx, double ny,
                       double ds, PRIM_VAR vari, PRIM_VAR varj, double gammai,
                       double gammaj, CONS_VAR &resi, CONS_VAR &resj) override;
};

class NumericSLAU2 : public BaseNumeric {
 public:
  using BaseNumeric::BaseNumeric;
  void DissipNumeric(const double &beta) override;
  void FluxNumeric() override;
  void ComputeResidual(const preprocess::parameter &param, double nx, double ny,
                       double ds, PRIM_VAR vari, PRIM_VAR varj, double gammai,
                       double gammaj, CONS_VAR &resi, CONS_VAR &resj) override;
};

class NumericAUSM : public BaseNumeric {
 public:
  using BaseNumeric::BaseNumeric;
  void DissipNumeric(const double &beta) override;
  void FluxNumeric() override;
  void ComputeResidual(const preprocess::parameter &param, double nx, double ny,
                       double ds, PRIM_VAR vari, PRIM_VAR varj, double gammai,
                       double gammaj, CONS_VAR &resi, CONS_VAR &resj) override;
};

class NumericAUSMUP2 : public BaseNumeric {
 public:
  using BaseNumeric::BaseNumeric;
  void DissipNumeric(const double &beta) override;
  void FluxNumeric() override;
  void ComputeResidual(const preprocess::parameter &param, double nx, double ny,
                       double ds, PRIM_VAR vari, PRIM_VAR varj, double gammai,
                       double gammaj, CONS_VAR &resi, CONS_VAR &resj) override;
};
}  // namespace solver