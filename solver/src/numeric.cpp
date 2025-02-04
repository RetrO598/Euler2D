#include "pre/parameter.h"
#include "solver/variableDef.h"
#include <algorithm>
#include <iostream>
#include <solver/numeric.h>

#include <cmath>

namespace solver {
BaseNumeric::BaseNumeric(const preprocess::parameter &param,
                         const preprocess::Geometry &geom,
                         std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
                         std::vector<CONS_VAR> &diss,
                         std::vector<CONS_VAR> &rhs, std::vector<PRIM_VAR> &lim,
                         std::vector<PRIM_VAR> &gradx,
                         std::vector<PRIM_VAR> &grady)
    : param(param), geom(geom), cv(cv), dv(dv), diss(diss), rhs(rhs), lim(lim),
      gradx(gradx), grady(grady) {}

void BaseNumeric::DissipInit(const int &irk, const double &beta) {
  if (irk == 0 || beta > 0.99) {
    for (int i = 0; i < geom.totNodes; ++i) {
      diss[i].dens = 0.0;
      diss[i].xmom = 0.0;
      diss[i].ymom = 0.0;
      diss[i].ener = 0.0;
    }
  } else {
    double blend = 1.0 - beta;
    for (int i = 0; i < geom.totNodes; ++i) {
      diss[i].dens *= blend;
      diss[i].xmom *= blend;
      diss[i].ymom *= blend;
      diss[i].ener *= blend;
    }
  }
}

void BaseNumeric::FluxWalls() {
  int ibegf = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendf = geom.ibound[ib].bfaceIndex;
    if (geom.BoundTypes[ib] >= 300 && geom.BoundTypes[ib] < 500) {
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

double NumericRoe::EntropyCorr(const double &z, const double &d) {
  if (z > d) {
    return z;
  } else {
    return 0.5 * (z * z + d * d) / d;
  }
}

void NumericRoe::DissipNumeric(const double &beta) {
  double beta5 = 0.5 * beta;
  double rrho, gam1, ggm1, rl, ul, vl, pl, hl, rr, ur, vr, pr, hr, rav, dd, dd1,
      uav, vav, hav, q2a, c2a, cav, uv, du, h1, h2, h3, h4, h5, delta, eabs1,
      eabs2, eabs4;
  double fd[4];

  for (int ie = 0; ie < geom.totEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double nx = geom.sij[ie].x / ds;
    double ny = geom.sij[ie].y / ds;

    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);

    interfaceVar vars = Interpolate(i, j);

    rl = vars.left.dens;
    ul = vars.left.velx;
    vl = vars.left.vely;
    pl = vars.left.press;
    hl =
        dv[i].gamma / (dv[i].gamma - 1.0) * pl / rl + 0.5 * (ul * ul + vl * vl);

    rr = vars.right.dens;
    ur = vars.right.velx;
    vr = vars.right.vely;
    pr = vars.right.press;
    hr =
        dv[j].gamma / (dv[j].gamma - 1.0) * pr / rr + 0.5 * (ur * ur + vr * vr);

    rav = std::sqrt(rl * rr);
    gam1 = 0.5 * (dv[i].gamma + dv[j].gamma) - 1.0;
    dd = rav / rl;
    dd1 = 1.0 / (1.0 + dd);
    uav = (ul + dd * ur) * dd1;
    vav = (vl + dd * vr) * dd1;
    hav = (hl + dd * hr) * dd1;
    q2a = 0.5 * (uav * uav + vav * vav);
    c2a = gam1 * (hav - q2a);
    cav = std::sqrt(c2a);
    uv = uav * nx + vav * ny;
    du = (ur - ul) * nx + (vr - vl) * ny;

    h1 = std::abs(uv - cav);
    h2 = std::abs(uv);
    h4 = std::abs(uv + cav);
    delta = param.entropyCorreRoe * h4;

    eabs1 = EntropyCorr(h1, delta);
    eabs2 = EntropyCorr(h2, delta);
    eabs4 = EntropyCorr(h4, delta);

    h1 = rav * cav * du;
    h2 = eabs1 * (pr - pl - h1) / (2.0 * c2a);
    h3 = eabs2 * (rr - rl - (pr - pl) / c2a);
    h4 = eabs2 * rav;
    h5 = eabs4 * (pr - pl + h1) / (2.0 * c2a);

    fd[0] = h2 + h3 + h5;
    fd[1] = h2 * (uav - cav * nx) + h3 * uav + h4 * (ur - ul - du * nx) +
            h5 * (uav + cav * nx);
    fd[2] = h2 * (vav - cav * ny) + h3 * vav + h4 * (vr - vl - du * ny) +
            h5 * (vav + cav * ny);
    fd[3] = h2 * (hav - cav * uv) + h3 * q2a +
            h4 * (uav * (ur - ul) + vav * (vr - vl) - uv * du) +
            h5 * (hav + cav * uv);

    ds *= beta5;
    diss[i].dens += fd[0] * ds;
    diss[i].xmom += fd[1] * ds;
    diss[i].ymom += fd[2] * ds;
    diss[i].ener += fd[3] * ds;

    diss[j].dens -= fd[0] * ds;
    diss[j].xmom -= fd[1] * ds;
    diss[j].ymom -= fd[2] * ds;
    diss[j].ener -= fd[3] * ds;
  }
}

void NumericRoe::FluxNumeric() {
  double rx, ry, rrho, gam1, ggm1, rl, ul, vl, pl, hl, rr, ur, vr, pr, hr, qsrl,
      qsrr, pav;
  double fc[4];

  for (int i = 0; i < geom.totNodes; ++i) {
    rhs[i].dens = -diss[i].dens;
    rhs[i].xmom = -diss[i].xmom;
    rhs[i].ymom = -diss[i].ymom;
    rhs[i].ener = -diss[i].ener;
  }

  for (int ie = 0; ie < geom.totEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);

    interfaceVar vars = Interpolate(i, j);

    rl = vars.left.dens;
    ul = vars.left.velx;
    vl = vars.left.vely;
    pl = vars.left.press;
    hl =
        dv[i].gamma / (dv[i].gamma - 1.0) * pl / rl + 0.5 * (ul * ul + vl * vl);

    rr = vars.right.dens;
    ur = vars.right.velx;
    vr = vars.right.vely;
    pr = vars.right.press;
    hr =
        dv[j].gamma / (dv[j].gamma - 1.0) * pr / rr + 0.5 * (ur * ur + vr * vr);

    qsrl = (ul * geom.sij[ie].x + vl * geom.sij[ie].y) * rl;
    qsrr = (ur * geom.sij[ie].x + vr * geom.sij[ie].y) * rr;

    pav = 0.5 * (pl + pr);
    fc[0] = 0.5 * (qsrl + qsrr);
    fc[1] = 0.5 * (qsrl * ul + qsrr * ur) + geom.sij[ie].x * pav;
    fc[2] = 0.5 * (qsrl * vl + qsrr * vr) + geom.sij[ie].y * pav;
    fc[3] = 0.5 * (qsrl * hl + qsrr * hr);

    rhs[i].dens += fc[0];
    rhs[i].xmom += fc[1];
    rhs[i].ymom += fc[2];
    rhs[i].ener += fc[3];

    rhs[j].dens -= fc[0];
    rhs[j].xmom -= fc[1];
    rhs[j].ymom -= fc[2];
    rhs[j].ener -= fc[3];
  }

  FluxWalls();
}

void NumericSLAU2::DissipNumeric(const double &beta) {
  for (int ie = 0; ie < geom.totEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    diss[i].dens += 0.0;
    diss[i].xmom += 0.0;
    diss[i].ymom += 0.0;
    diss[i].ener += 0.0;

    diss[j].dens -= 0.0;
    diss[j].xmom -= 0.0;
    diss[j].ymom -= 0.0;
    diss[j].ener -= 0.0;
  }
}

void NumericSLAU2::FluxNumeric() {
  double fc[4];
  for (int i = 0; i < geom.totNodes; ++i) {
    rhs[i].dens = -diss[i].dens;
    rhs[i].xmom = -diss[i].xmom;
    rhs[i].ymom = -diss[i].ymom;
    rhs[i].ener = -diss[i].ener;
  }
  for (int ie = 0; ie < geom.totEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double nx = geom.sij[ie].x / ds;
    double ny = geom.sij[ie].y / ds;
    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);

    interfaceVar vars = Interpolate(i, j);

    double rl = vars.left.dens;
    double ul = vars.left.velx;
    double vl = vars.left.vely;
    double pl = vars.left.press;
    double hl =
        dv[i].gamma / (dv[i].gamma - 1.0) * pl / rl + 0.5 * (ul * ul + vl * vl);

    double rr = vars.right.dens;
    double ur = vars.right.velx;
    double vr = vars.right.vely;
    double pr = vars.right.press;
    double hr =
        dv[j].gamma / (dv[j].gamma - 1.0) * pr / rr + 0.5 * (ur * ur + vr * vr);

    double cl = (dv[i].gamma - 1.0) * (hl - 0.5 * (ul * ul + vl * vl));
    double cr = (dv[j].gamma - 1.0) * (hr - 0.5 * (ur * ur + vr * vr));

    double veln_left = ul * nx + vl * ny;
    double veln_right = ur * nx + vr * ny;
    double veln_bar_abs =
        (rl * std::abs(veln_left) + rr * std::abs(veln_right)) / (rl + rr);
    double c_interface = 0.5 * (cl + cr);
    double m_left = veln_left / c_interface;
    double m_right = veln_right / c_interface;
    double g = -1.0 * std::max(std::min(m_left, 0.0), -1.0) *
               std::min(std::max(m_right, 0.0), 1.0);
    double veln_bar_abs_p = (1.0 - g) * veln_bar_abs + g * std::abs(veln_left);
    double veln_bar_abs_n = (1.0 - g) * veln_bar_abs + g * std::abs(veln_right);

    double m_hat = std::min(
        1.0, 1 / c_interface *
                 std::sqrt(0.5 * (ul * ul + vl * vl + ur * ur + vr * vr)));
    double X = std::pow((1.0 - m_hat), 2);

    double fp_left = pressureFuncLeft(m_left);
    double fp_right = pressureFuncRight(m_right);
    double p_interface =
        0.5 * (pl + pr) + 0.5 * (fp_left - fp_right) * (pl - pr) +
        std::sqrt(0.5 * (ul * ul + vl * vl + ur * ur + vr * vr)) *
            (fp_left + fp_right - 1) * c_interface * 0.5 * (rr + rl);
    double mass_interface = 0.5 * (rl * (veln_left + veln_bar_abs_p) +
                                   rr * (veln_right - veln_bar_abs_n) -
                                   X / c_interface * (pr - pl));
    double mass_flux_left = 0.5 * (mass_interface + std::abs(mass_interface));
    double mass_flux_right = 0.5 * (mass_interface - std::abs(mass_interface));

    fc[0] = (mass_flux_left + mass_flux_right) * ds;
    fc[1] =
        (mass_flux_left * ul + mass_flux_right * ur + p_interface * nx) * ds;
    fc[2] =
        (mass_flux_left * vl + mass_flux_right * vr + p_interface * ny) * ds;
    fc[3] = (mass_flux_left * hl + mass_flux_right * hr) * ds;

    rhs[i].dens += fc[0];
    rhs[i].xmom += fc[1];
    rhs[i].ymom += fc[2];
    rhs[i].ener += fc[3];

    rhs[j].dens -= fc[0];
    rhs[j].xmom -= fc[1];
    rhs[j].ymom -= fc[2];
    rhs[j].ener -= fc[3];
  }
  FluxWalls();
}

void NumericAUSM::DissipNumeric(const double &beta) {
  for (int ie = 0; ie < geom.totEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    diss[i].dens += 0.0;
    diss[i].xmom += 0.0;
    diss[i].ymom += 0.0;
    diss[i].ener += 0.0;

    diss[j].dens -= 0.0;
    diss[j].xmom -= 0.0;
    diss[j].ymom -= 0.0;
    diss[j].ener -= 0.0;
  }
}

void NumericAUSM::FluxNumeric() {
  double fc[4];
  for (int i = 0; i < geom.totNodes; ++i) {
    rhs[i].dens = -diss[i].dens;
    rhs[i].xmom = -diss[i].xmom;
    rhs[i].ymom = -diss[i].ymom;
    rhs[i].ener = -diss[i].ener;
  }
  for (int ie = 0; ie < geom.totEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double nx = geom.sij[ie].x / ds;
    double ny = geom.sij[ie].y / ds;
    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);

    interfaceVar vars = Interpolate(i, j);

    double rl = vars.left.dens;
    double ul = vars.left.velx;
    double vl = vars.left.vely;
    double pl = vars.left.press;
    double hl =
        dv[i].gamma / (dv[i].gamma - 1.0) * pl / rl + 0.5 * (ul * ul + vl * vl);

    double rr = vars.right.dens;
    double ur = vars.right.velx;
    double vr = vars.right.vely;
    double pr = vars.right.press;
    double hr =
        dv[j].gamma / (dv[j].gamma - 1.0) * pr / rr + 0.5 * (ur * ur + vr * vr);
    double cl2 = (dv[i].gamma - 1.0) * (hl - 0.5 * (ul * ul + vl * vl));
    double cl = std::sqrt(cl2);
    double cr2 = (dv[j].gamma - 1.0) * (hr - 0.5 * (ur * ur + vr * vr));
    double cr = std::sqrt(cr2);

    double veln_left = ul * nx + vl * ny;
    double veln_right = ur * nx + vr * ny;

    double m_left = veln_left / cl;
    double m_right = veln_right / cr;

    double m_left_p;
    double p_left_p;
    if (m_left >= 1.0) {
      m_left_p = m_left;
      p_left_p = pl;
    } else if (std::abs(m_left) < 1.0) {
      m_left_p = 0.25 * (m_left + 1.0) * (m_left + 1.0);
      p_left_p = 0.25 * pl * (m_left + 1.0) * (m_left + 1.0) * (2.0 - m_left);
    } else if (m_left <= -1.0) {
      m_left_p = 0.0;
      p_left_p = 0.0;
    }

    double m_right_n;
    double p_right_n;
    if (m_right >= 1.0) {
      m_right_n = 0.0;
      p_right_n = 0.0;
    } else if (std::abs(m_right) < 1.0) {
      m_right_n = -0.25 * (m_right - 1.0) * (m_right - 1.0);
      p_right_n =
          0.25 * pr * (m_right - 1.0) * (m_right - 1.0) * (2.0 + m_right);
    } else if (m_right <= -1.0) {
      m_right_n = m_right;
      p_right_n = pr;
    }

    double m_interface = m_left_p + m_right_n;
    double p_interface = p_left_p + p_right_n;
    double m_interface_abs;
    double delta = 0.1;
    if (std::abs(m_interface) > delta) {
      m_interface_abs = std::abs(m_interface);
    } else {
      m_interface_abs =
          (m_interface * m_interface + delta * delta) / (2.0 * delta);
    }
    fc[0] = 0.5 * m_interface * (rl * cl + rr * cr) * ds -
            0.5 * m_interface_abs * (rr * cr - rl * cl) * ds;
    fc[1] = 0.5 * m_interface * (rl * cl * ul + rr * cr * ur) * ds -
            0.5 * m_interface_abs * (rr * cr * ur - rl * cl * ul) * ds +
            p_interface * nx * ds;
    fc[2] = 0.5 * m_interface * (rl * cl * vl + rr * cr * vr) * ds -
            0.5 * m_interface_abs * (rr * cr * vr - rl * cl * vl) * ds +
            p_interface * ny * ds;
    fc[3] = 0.5 * m_interface * (rl * cl * hl + rr * cr * hr) * ds -
            0.5 * m_interface_abs * (rr * cr * hr - rl * cl * hl) * ds;

    rhs[i].dens += fc[0];
    rhs[i].xmom += fc[1];
    rhs[i].ymom += fc[2];
    rhs[i].ener += fc[3];

    rhs[j].dens -= fc[0];
    rhs[j].xmom -= fc[1];
    rhs[j].ymom -= fc[2];
    rhs[j].ener -= fc[3];
  }
  FluxWalls();
}

void NumericAUSMUP2::DissipNumeric(const double &beta) {
  for (int ie = 0; ie < geom.totEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    diss[i].dens += 0.0;
    diss[i].xmom += 0.0;
    diss[i].ymom += 0.0;
    diss[i].ener += 0.0;

    diss[j].dens -= 0.0;
    diss[j].xmom -= 0.0;
    diss[j].ymom -= 0.0;
    diss[j].ener -= 0.0;
  }
}

void NumericAUSMUP2::FluxNumeric() {
  double fc[4];
  for (int i = 0; i < geom.totNodes; ++i) {
    rhs[i].dens = -diss[i].dens;
    rhs[i].xmom = -diss[i].xmom;
    rhs[i].ymom = -diss[i].ymom;
    rhs[i].ener = -diss[i].ener;
  }
  for (int ie = 0; ie < geom.totEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;
    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double nx = geom.sij[ie].x / ds;
    double ny = geom.sij[ie].y / ds;
    double rx = 0.5 * (geom.coords[j].x - geom.coords[i].x);
    double ry = 0.5 * (geom.coords[j].y - geom.coords[i].y);

    interfaceVar vars = Interpolate(i, j);

    double rl = vars.left.dens;
    double ul = vars.left.velx;
    double vl = vars.left.vely;
    double pl = vars.left.press;
    double hl =
        dv[i].gamma / (dv[i].gamma - 1.0) * pl / rl + 0.5 * (ul * ul + vl * vl);

    double rr = vars.right.dens;
    double ur = vars.right.velx;
    double vr = vars.right.vely;
    double pr = vars.right.press;
    double hr =
        dv[j].gamma / (dv[j].gamma - 1.0) * pr / rr + 0.5 * (ur * ur + vr * vr);
    double cl2 = 2.0 * (dv[i].gamma - 1.0) / (dv[i].gamma + 1.0) * hl;
    double cl = std::sqrt(cl2);
    double cr2 = 2.0 * (dv[j].gamma - 1.0) / (dv[j].gamma + 1.0) * hr;
    double cr = std::sqrt(cr2);

    double veln_left = ul * nx + vl * ny;
    double veln_right = ur * nx + vr * ny;
    double c_bar_left = cl2 / (std::max(cl, veln_left));
    double c_bar_right = cr2 / (std::max(cr, -veln_right));

    double c_interface = std::min(c_bar_left, c_bar_right);
    double m_left = veln_left / c_interface;
    double m_right = veln_right / c_interface;
    double m_bar_2 = (veln_left * veln_left + veln_right * veln_right) /
                     (2.0 * c_interface * c_interface);
    double m_infty;
    if (param.flowtype_ == preprocess::flowType::External) {
      m_infty = param.MaInfinity;
    } else {
      double mach = std::sqrt(param.refMach2);
      m_infty = 0.2 * std::min(mach, 1.0);
    }
    double m0_2 = std::min(1.0, std::max(m_bar_2, m_infty * m_infty));
    double m0 = std::sqrt(m0_2);
    double fa = m0 * (2.0 - m0);
    double mp = -0.25 / fa * std::max(1.0 - m_bar_2, 0.0) * (pr - pl) /
                (0.5 * (rr + rl) * c_interface * c_interface);

    double fp_left;
    double fp_right;
    double fm_left;
    double fm_right;
    if (std::abs(m_left) >= 1.0) {
      fm_left = 0.5 * (m_left + std::abs(m_left));
      fp_left = 0.5 * (1.0 + (m_left >= 0.0 ? 1.0 : -1.0));
    } else {
      fm_left = 0.25 * (m_left + 1.0) * (m_left + 1.0) +
                0.125 * (m_left * m_left - 1.0) * (m_left * m_left - 1.0);
      fp_left = 0.25 * (m_left + 1.0) * (m_left + 1.0) * (2.0 - m_left);
    }
    if (std::abs(m_right) >= 1.0) {
      fm_right = 0.5 * (m_right - std::abs(m_right));
      fp_right = 0.5 * (1.0 - (m_right >= 0.0 ? 1.0 : -1.0));
    } else {
      fm_right = -0.25 * (m_right - 1.0) * (m_right - 1.0) -
                 0.125 * (m_right * m_right - 1.0) * (m_right * m_right - 1.0);
      fp_right = 0.25 * (m_right - 1.0) * (m_right - 1.0) * (2.0 + m_right);
    }
    double m_interface = fm_left + fm_right + mp;

    double p_interface =
        0.5 * (pl + pr) + 0.5 * (fp_left - fp_right) * (pl - pr) +
        std::sqrt(0.5 * (ul * ul + vl * vl + ur * ur + vr * vr)) *
            (fp_left + fp_right - 1) * c_interface * 0.5 * (rr + rl);

    double mass_dot;
    if (m_interface > 0.0) {
      mass_dot = m_interface * c_interface * rl;
    } else {
      mass_dot = m_interface * c_interface * rr;
    }

    double mass_flux_left = 0.5 * (mass_dot + std::abs(mass_dot));
    double mass_flux_right = 0.5 * (mass_dot - std::abs(mass_dot));

    fc[0] = (mass_flux_left + mass_flux_right) * ds;
    fc[1] =
        (mass_flux_left * ul + mass_flux_right * ur + p_interface * nx) * ds;
    fc[2] =
        (mass_flux_left * vl + mass_flux_right * vr + p_interface * ny) * ds;
    fc[3] = (mass_flux_left * hl + mass_flux_right * hr) * ds;

    rhs[i].dens += fc[0];
    rhs[i].xmom += fc[1];
    rhs[i].ymom += fc[2];
    rhs[i].ener += fc[3];

    rhs[j].dens -= fc[0];
    rhs[j].xmom -= fc[1];
    rhs[j].ymom -= fc[2];
    rhs[j].ener -= fc[3];
  }
  FluxWalls();
}
} // namespace solver