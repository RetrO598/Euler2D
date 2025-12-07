#pragma once

#include <sys/stat.h>

#include <functional>
#include <memory>
#include <unordered_map>

#include "pre/parameter.h"
#include "solver/limiter.h"
#include "solver/numeric.h"
namespace solver {
using NumericPtr = std::unique_ptr<BaseNumeric>;
using NumericCreatorFunc = std::function<NumericPtr(
    const preprocess::parameter &param, const preprocess::Geometry &geom,
    std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
    std::vector<CONS_VAR> &diss, std::vector<CONS_VAR> &rhs,
    std::vector<PRIM_VAR> &lim, std::vector<PRIM_VAR> &gradx,
    std::vector<PRIM_VAR> &grady)>;

class ConvectionFactory {
 public:
  static NumericPtr create(
      const preprocess::parameter &param, const preprocess::Geometry &geom,
      std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
      std::vector<CONS_VAR> &diss, std::vector<CONS_VAR> &rhs,
      std::vector<PRIM_VAR> &lim, std::vector<PRIM_VAR> &gradx,
      std::vector<PRIM_VAR> &grady);

 private:
  static std::unordered_map<preprocess::ConvectionScheme, NumericCreatorFunc> &
  getRegistry();
};

using LimiterPtr = std::unique_ptr<BaseLimiter>;
using LimiterCreatorFunc = std::function<LimiterPtr(
    const preprocess::parameter &param, const preprocess::Geometry &geom,
    std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
    std::vector<PRIM_VAR> &umin, std::vector<PRIM_VAR> &umax,
    std::vector<PRIM_VAR> &lim, std::vector<PRIM_VAR> &gradx,
    std::vector<PRIM_VAR> &grady)>;

class LimiterFactory {
 public:
  static LimiterPtr create(
      const preprocess::parameter &param, const preprocess::Geometry &geom,
      std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
      std::vector<PRIM_VAR> &umin, std::vector<PRIM_VAR> &umax,
      std::vector<PRIM_VAR> &lim, std::vector<PRIM_VAR> &gradx,
      std::vector<PRIM_VAR> &grady);

 private:
  static std::unordered_map<preprocess::Limiter, LimiterCreatorFunc> &
  getRegistry();
};
}  // namespace solver