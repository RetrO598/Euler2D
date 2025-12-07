#include <iostream>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "pre/parameter.h"
#include "solver/limiter.h"
#include "solver/numeric.h"
#include "solver/numericFactory.hpp"
namespace solver {
NumericPtr ConvectionFactory::create(
    const preprocess::parameter &param, const preprocess::Geometry &geom,
    std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
    std::vector<CONS_VAR> &diss, std::vector<CONS_VAR> &rhs,
    std::vector<PRIM_VAR> &lim, std::vector<PRIM_VAR> &gradx,
    std::vector<PRIM_VAR> &grady) {
  const auto &registry = getRegistry();
  auto it = registry.find(param.convecScheme);
  if (it == registry.end()) {
    throw std::runtime_error("unknown convection scheme.\n");
  }
  return it->second(param, geom, cv, dv, diss, rhs, lim, gradx, grady);
}

std::unordered_map<preprocess::ConvectionScheme, NumericCreatorFunc> &
ConvectionFactory::getRegistry() {
  static std::unordered_map<preprocess::ConvectionScheme, NumericCreatorFunc>
      registry = {{preprocess::ConvectionScheme::ROE,
                   [](auto &&...args) {
                     std::cout << "using ROE scheme.\n";
                     return std::make_unique<NumericRoe>(
                         std::forward<decltype(args)>(args)...);
                   }},
                  {preprocess::ConvectionScheme::SLAU2,
                   [](auto &&...args) {
                     std::cout << "using SLAU2 scheme.\n";
                     return std::make_unique<NumericSLAU2>(
                         std::forward<decltype(args)>(args)...);
                   }},
                  {preprocess::ConvectionScheme::AUSM,
                   [](auto &&...args) {
                     std::cout << "using AUSM scheme.\n";
                     return std::make_unique<NumericAUSM>(
                         std::forward<decltype(args)>(args)...);
                   }},
                  {preprocess::ConvectionScheme::AUSMUP2, [](auto &&...args) {
                     std::cout << "using AUSMUP2 scheme.\n";
                     return std::make_unique<NumericAUSMUP2>(
                         std::forward<decltype(args)>(args)...);
                   }}};
  return registry;
}

LimiterPtr LimiterFactory::create(
    const preprocess::parameter &param, const preprocess::Geometry &geom,
    std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
    std::vector<PRIM_VAR> &umin, std::vector<PRIM_VAR> &umax,
    std::vector<PRIM_VAR> &lim, std::vector<PRIM_VAR> &gradx,
    std::vector<PRIM_VAR> &grady) {
  const auto &registry = getRegistry();
  auto it = registry.find(param.limiterType);
  if (it == registry.end()) {
    throw std::runtime_error("unknown limiter type.\n");
  }
  return it->second(param, geom, cv, dv, umin, umax, lim, gradx, grady);
}

std::unordered_map<preprocess::Limiter, LimiterCreatorFunc> &
LimiterFactory::getRegistry() {
  static std::unordered_map<preprocess::Limiter, LimiterCreatorFunc> registry =
      {{preprocess::Limiter::VenkataKrishnan,
        [](auto &&...args) {
          std::cout << "using Venkatakrishnan limiter.\n";
          return std::make_unique<VenkatakrishnanLimiter>(
              std::forward<decltype(args)>(args)...);
        }},
       {preprocess::Limiter::NishikawaR3, [](auto &&...args) {
          std::cout << "using NishikawaR3 limiter.\n";
          return std::make_unique<NishikawaR3>(
              std::forward<decltype(args)>(args)...);
        }}};
  return registry;
}
}  // namespace solver