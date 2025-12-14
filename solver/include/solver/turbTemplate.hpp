#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "solver/variableDef.h"
namespace solver {
struct UpwindTurbSA {
  inline static void ComputeResidual(double nx, double ny, double ds,
                                     TurbSA_VAR nui, TurbSA_VAR nuj,
                                     CONS_VAR vari, CONS_VAR varj,
                                     TurbSA_VAR &resi, TurbSA_VAR &resj) {
    double ul = vari.xmom / vari.dens;
    double vl = vari.ymom / vari.dens;
    double ur = varj.xmom / varj.dens;
    double vr = varj.ymom / varj.dens;

    double q = 0.5 * (ul + ur) * nx + 0.5 * (vl + vr) * ny;
    double a0 = 0.5 * (q + std::fabs(q));
    double a1 = 0.5 * (q - std::fabs(q));

    double result = a0 * nui.nu_turb + a1 * nuj.nu_turb;

    resi.nu_turb += result * ds;
    resj.nu_turb -= result * ds;

    if (std::isnan(result)) {
      std::cout << "not a number from convection\n";
      // exit(1);
    }
  }
};

struct TurbViscousBound {
  inline static void ComputeResidual(double nuLami, double nuLamj,
                                     double nuTurbi, double nuTurbj,
                                     double gradx, double grady, double sx,
                                     double sy, TurbSA_VAR &resi) {
    double nuLamAve = 0.5 * (nuLami + nuLamj);
    double nuTurbAve = 0.5 * (nuTurbi + nuTurbj);
    double sigma = 2.0 / 3.0;
    double fv;
    fv = (nuLamAve + nuTurbAve) / sigma * gradx * sx +
         (nuLamAve + nuTurbAve) / sigma * grady * sy;

    resi.nu_turb -= fv;

    if (std::isnan(fv)) {
      std::cout << "not a number from visous boundary\n";
      // exit(1);
    }
  }
};
}  // namespace solver