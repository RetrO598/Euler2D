#pragma once
#include <cmath>
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
  void limiterRefVals();
  void limiterInit();
  virtual void limiterUpdate(){};

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
  using BaseLimiter::BaseLimiter;
  void limiterUpdate() override;

  inline double Venkat(double d2, double d1min, double d1max, double eps2) {
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

class NishikawaR3 : public BaseLimiter {
public:
  using BaseLimiter::BaseLimiter;
  void limiterUpdate() override;

  inline double Nishikawa_R3(const double &delta_pos, const double &delta_neg,
                             const double &eps) {
    double a = std::abs(delta_pos);
    double b = std::abs(delta_neg);

    if (a > 2.0 * b) {
      return 1.0;
    } else {
      double ap = std::pow(a, 3);
      double bp = std::pow(b, 3);
      double Sp = 4.0 * b * b;
      double num = ap + eps + a * Sp;
      double den = ap + eps + b * (std::pow(delta_pos, 2) + Sp);
      return num / den;
    }
  }
};
} // namespace solver