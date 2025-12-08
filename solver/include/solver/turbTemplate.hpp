#pragma once

#include <cmath>

#include "solver/variableDef.h"
namespace solver {
struct UpwindTurbSA {
  inline static void ComputeResidual(double nx, double ny, double ds,
                                     TurbSA_VAR nui, TurbSA_VAR nuj,
                                     PRIM_VAR vari, PRIM_VAR varj,
                                     TurbSA_VAR &resi, TurbSA_VAR &resj) {
    double ul = vari.velx;
    double vl = vari.vely;
    double ur = varj.velx;
    double vr = varj.vely;

    double q = 0.5 * (ul + ur) * nx + 0.5 * (vl + vr) * ny;
    double a0 = 0.5 * (q + std::fabs(q));
    double a1 = 0.5 * (q - std::fabs(q));

    double result = a0 * nui.nu_turb + a1 * nuj.nu_turb;

    resi.nu_turb += result * ds;
    resj.nu_turb -= result * ds;
  }
};
}  // namespace solver