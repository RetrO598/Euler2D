#include <solver/numeric.h>

#include <cmath>
#include <solver/numericTemplate.hpp>

#include "pre/parameter.h"
#include "solver/variableDef.h"

namespace solver {
BaseNumeric::BaseNumeric(const preprocess::parameter &param,
                         const preprocess::Geometry &geom,
                         std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
                         std::vector<CONS_VAR> &diss,
                         std::vector<CONS_VAR> &rhs, std::vector<PRIM_VAR> &lim,
                         std::vector<PRIM_VAR> &gradx,
                         std::vector<PRIM_VAR> &grady)
    : param(param),
      geom(geom),
      cv(cv),
      dv(dv),
      diss(diss),
      rhs(rhs),
      lim(lim),
      gradx(gradx),
      grady(grady) {}

// void BaseNumeric::DissipInit(const int &irk, const double &beta) {
//   if (irk == 0 || beta > 0.99) {
//     for (int i = 0; i < geom.phyNodes; ++i) {
//       diss[i].dens = 0.0;
//       diss[i].xmom = 0.0;
//       diss[i].ymom = 0.0;
//       diss[i].ener = 0.0;
//     }
//   } else {
//     double blend = 1.0 - beta;
//     for (int i = 0; i < geom.phyNodes; ++i) {
//       diss[i].dens *= blend;
//       diss[i].xmom *= blend;
//       diss[i].ymom *= blend;
//       diss[i].ener *= blend;
//     }
//   }
// }

void BaseNumeric::DissipInit() {
  for (int i = 0; i < geom.phyNodes; ++i) {
    diss[i].dens = 0.0;
    diss[i].xmom = 0.0;
    diss[i].ymom = 0.0;
    diss[i].ener = 0.0;
  }
}

void BaseNumeric::FluxWalls() {
  int ibegf = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendf = geom.ibound[ib].bfaceIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::EulerWall ||
        type == preprocess::BoundaryType::NoSlipWall) {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        int i = geom.boundaryFace[ibf].nodei;
        int j = geom.boundaryFace[ibf].nodej;
        double sx = geom.sbf[ibf].x / 12.0;
        double sy = geom.sbf[ibf].y / 12.0;
        double pl = 5.0 * dv[i].press + dv[j].press;
        double pr = dv[i].press + 5.0 * dv[j].press;
        rhs[i].xmom += sx * pl;
        rhs[i].ymom += sy * pl;
        rhs[j].xmom += sx * pr;
        rhs[j].ymom += sy * pr;
      }
    }
    ibegf = iendf + 1;
  }
}

interfaceVar BaseNumeric::Interpolate(const int &i, const int &j) {
  double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
  double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);
  double rrho = 1.0 / cv[i].dens;
  double gam1 = dv[i].gamma - 1.0;
  double ggm1 = dv[i].gamma / gam1;
  double r =
      cv[i].dens + lim[i].dens * (gradx[i].dens * rx + grady[i].dens * ry);
  double u = cv[i].xmom * rrho +
             lim[i].velx * (gradx[i].velx * rx + grady[i].velx * ry);
  double v = cv[i].ymom * rrho +
             lim[i].vely * (gradx[i].vely * rx + grady[i].vely * ry);
  double p =
      dv[i].press + lim[i].press * (gradx[i].press * rx + grady[i].press * ry);

  PRIM_VAR left{r, u, v, p};

  rrho = 1.0 / cv[j].dens;
  gam1 = dv[j].gamma - 1.0;
  ggm1 = dv[j].gamma / gam1;
  r = cv[j].dens - lim[j].dens * (gradx[j].dens * rx + grady[j].dens * ry);
  u = cv[j].xmom * rrho -
      lim[j].velx * (gradx[j].velx * rx + grady[j].velx * ry);
  v = cv[j].ymom * rrho -
      lim[j].vely * (gradx[j].vely * rx + grady[j].vely * ry);
  p = dv[j].press - lim[j].press * (gradx[j].press * rx + grady[j].press * ry);
  PRIM_VAR right{r, u, v, p};

  return interfaceVar{left, right};
}

void NumericRoe::ComputeResidual(const preprocess::parameter &param, double nx,
                                 double ny, double ds, PRIM_VAR vari,
                                 PRIM_VAR varj, double gammai, double gammaj,
                                 CONS_VAR &resi, CONS_VAR &resj) {
  ROENumeric::ComputeResidual(param, nx, ny, ds, vari, varj, gammai, gammaj,
                              resi, resj);
};

void NumericRoe::FluxNumeric() {
  for (int i = 0; i < geom.phyNodes; ++i) {
    rhs[i].dens = -diss[i].dens;
    rhs[i].xmom = -diss[i].xmom;
    rhs[i].ymom = -diss[i].ymom;
    rhs[i].ener = -diss[i].ener;
  }

  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double nx = geom.sij[ie].x / ds;
    double ny = geom.sij[ie].y / ds;

    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);

    interfaceVar vars = Interpolate(i, j);

    ComputeResidual(param, nx, ny, ds, vars.left, vars.right, dv[i].gamma,
                    dv[j].gamma, rhs[i], rhs[j]);
  }

  FluxWalls();
}

void NumericSLAU2::ComputeResidual(const preprocess::parameter &param,
                                   double nx, double ny, double ds,
                                   PRIM_VAR vari, PRIM_VAR varj, double gammai,
                                   double gammaj, CONS_VAR &resi,
                                   CONS_VAR &resj) {
  SLAU2Numeric::ComputeResidual(param, nx, ny, ds, vari, varj, gammai, gammaj,
                                resi, resj);
}

void NumericSLAU2::FluxNumeric() {
  double fc[4];
  for (int i = 0; i < geom.phyNodes; ++i) {
    rhs[i].dens = -diss[i].dens;
    rhs[i].xmom = -diss[i].xmom;
    rhs[i].ymom = -diss[i].ymom;
    rhs[i].ener = -diss[i].ener;
  }
  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double nx = geom.sij[ie].x / ds;
    double ny = geom.sij[ie].y / ds;
    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);

    interfaceVar vars = Interpolate(i, j);

    ComputeResidual(param, nx, ny, ds, vars.left, vars.right, dv[i].gamma,
                    dv[j].gamma, rhs[i], rhs[j]);
  }
  FluxWalls();
}

void NumericAUSM::ComputeResidual(const preprocess::parameter &param, double nx,
                                  double ny, double ds, PRIM_VAR vari,
                                  PRIM_VAR varj, double gammai, double gammaj,
                                  CONS_VAR &resi, CONS_VAR &resj) {
  AUSMNumeric::ComputeResidual(param, nx, ny, ds, vari, varj, gammai, gammaj,
                               resi, resj);
}

void NumericAUSM::FluxNumeric() {
  double fc[4];
  for (int i = 0; i < geom.phyNodes; ++i) {
    rhs[i].dens = -diss[i].dens;
    rhs[i].xmom = -diss[i].xmom;
    rhs[i].ymom = -diss[i].ymom;
    rhs[i].ener = -diss[i].ener;
  }
  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double nx = geom.sij[ie].x / ds;
    double ny = geom.sij[ie].y / ds;
    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);

    interfaceVar vars = Interpolate(i, j);

    ComputeResidual(param, nx, ny, ds, vars.left, vars.right, dv[i].gamma,
                    dv[j].gamma, rhs[i], rhs[j]);
  }
  FluxWalls();
}

void NumericAUSMUP2::ComputeResidual(const preprocess::parameter &param,
                                     double nx, double ny, double ds,
                                     PRIM_VAR vari, PRIM_VAR varj,
                                     double gammai, double gammaj,
                                     CONS_VAR &resi, CONS_VAR &resj) {
  if (param.flowtype_ == preprocess::flowType::External) {
    AUSMUP2Numeric::ComputeResidual<preprocess::flowType::External>(
        param, nx, ny, ds, vari, varj, gammai, gammaj, resi, resj);
  } else if (param.flowtype_ == preprocess::flowType::Internal) {
    AUSMUP2Numeric::ComputeResidual<preprocess::flowType::Internal>(
        param, nx, ny, ds, vari, varj, gammai, gammaj, resi, resj);
  }
}

void NumericAUSMUP2::FluxNumeric() {
  double fc[4];
  for (int i = 0; i < geom.phyNodes; ++i) {
    rhs[i].dens = -diss[i].dens;
    rhs[i].xmom = -diss[i].xmom;
    rhs[i].ymom = -diss[i].ymom;
    rhs[i].ener = -diss[i].ener;
  }
  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double nx = geom.sij[ie].x / ds;
    double ny = geom.sij[ie].y / ds;
    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);

    interfaceVar vars = Interpolate(i, j);

    ComputeResidual(param, nx, ny, ds, vars.left, vars.right, dv[i].gamma,
                    dv[j].gamma, rhs[i], rhs[j]);
  }
  FluxWalls();
}
}  // namespace solver