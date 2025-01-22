#pragma once
namespace solver {
struct CONS_VAR {
  double dens;
  double xmom;
  double ymom;
  double ener;
};

struct PRIM_VAR {
  double dens;
  double velx;
  double vely;
  double press;
};

struct DEPEND_VAR {
  double press;
  double temp;
  double cs;
  double gamma;
  double Cp;
};

struct VISC_VAR {
  double mu;
  double lambda;
};
}  // namespace solver