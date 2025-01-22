#pragma once
#include <pre/geometry.h>
#include <pre/parameter.h>
#include <solver/variableDef.h>

#include <vector>

namespace solver {

class limiter {
 public:
  virtual void limiterRefVals(preprocess::parameter &param,
                              const preprocess::Geometry &geom) = 0;
  virtual void limiterInit(preprocess::parameter &param,
                           const preprocess::Geometry &geom,
                           std::vector<CONS_VAR> &cv,
                           std::vector<DEPEND_VAR> &dv,
                           std::vector<PRIM_VAR> &umin,
                           std::vector<PRIM_VAR> &umax,
                           std::vector<PRIM_VAR> &lim) = 0;
  virtual void limiterUpdate(
      preprocess::parameter &param, const preprocess::Geometry &geom,
      std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
      std::vector<PRIM_VAR> &umin, std::vector<PRIM_VAR> &umax,
      std::vector<PRIM_VAR> &lim, std::vector<PRIM_VAR> &gradx,
      std::vector<PRIM_VAR> &grady) = 0;
};

class VenkatakrishnanLimiter : limiter {
 public:
  VenkatakrishnanLimiter();
  void limiterRefVals(preprocess::parameter &param,
                      const preprocess::Geometry &geom) override;
  void limiterInit(preprocess::parameter &param,
                   const preprocess::Geometry &geom, std::vector<CONS_VAR> &cv,
                   std::vector<DEPEND_VAR> &dv, std::vector<PRIM_VAR> &umin,
                   std::vector<PRIM_VAR> &umax,
                   std::vector<PRIM_VAR> &lim) override;
  void limiterUpdate(preprocess::parameter &param,
                     const preprocess::Geometry &geom,
                     std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
                     std::vector<PRIM_VAR> &umin, std::vector<PRIM_VAR> &umax,
                     std::vector<PRIM_VAR> &lim, std::vector<PRIM_VAR> &gradx,
                     std::vector<PRIM_VAR> &grady) override;

  double Venkat(double d2, double d1min, double d1max, double eps2);

 public:
  PRIM_VAR limRef;
  double volRef;
};
}  // namespace solver