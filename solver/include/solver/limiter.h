#pragma once
#include <pre/geometry.h>
#include <pre/parameter.h>
#include <solver/variableDef.h>

#include <vector>

namespace solver {

class BaseLimiter {
 public:
  BaseLimiter(const preprocess::parameter &param,
              const preprocess::Geometry &geom, std::vector<CONS_VAR> &cv,
              std::vector<DEPEND_VAR> &dv, std::vector<PRIM_VAR> &umin,
              std::vector<PRIM_VAR> &umax, std::vector<PRIM_VAR> &lim,
              std::vector<PRIM_VAR> &gradx, std::vector<PRIM_VAR> &grady);
  virtual void limiterRefVals() {};
  virtual void limiterInit() {};
  virtual void limiterUpdate() {};

 protected:
  std::vector<CONS_VAR> &cv;
  std::vector<DEPEND_VAR> &dv;
  std::vector<PRIM_VAR> &umin;
  std::vector<PRIM_VAR> &umax;
  std::vector<PRIM_VAR> &lim;
  std::vector<PRIM_VAR> &gradx;
  std::vector<PRIM_VAR> &grady;

  PRIM_VAR limRef;
  double volRef;

  const preprocess::parameter &param;
  const preprocess::Geometry &geom;
};

class VenkatakrishnanLimiter : public BaseLimiter {
 public:
  VenkatakrishnanLimiter(const preprocess::parameter &param,
                         const preprocess::Geometry &geom,
                         std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
                         std::vector<PRIM_VAR> &umin,
                         std::vector<PRIM_VAR> &umax,
                         std::vector<PRIM_VAR> &lim,
                         std::vector<PRIM_VAR> &gradx,
                         std::vector<PRIM_VAR> &grady);
  void limiterRefVals() override;
  void limiterInit() override;
  void limiterUpdate() override;

  double Venkat(double d2, double d1min, double d1max, double eps2) {
    if (d2 > 1e-16) {
      double num = (d1max * d1max + eps2) * d2 + 2.0 * d2 * d2 * d1max;
      double den = d2 * (d1max * d1max + 2.0 * d2 * d2 + d1max * d2 + eps2);
      return (num / den);
    } else if (d2 < -1e-16) {
      double num = (d1min * d1min + eps2) * d2 + 2.0 * d2 * d2 * d1min;
      double den = d2 * (d1min * d1min + 2.0 * d2 * d2 + d1min * d2 + eps2);
      return (num / den);
    }
    return 1.0;
  }
};
}  // namespace solver