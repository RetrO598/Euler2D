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

NumericRoe::NumericRoe(const preprocess::parameter &param,
                       const preprocess::Geometry &geom,
                       std::vector<CONS_VAR> &cv, std::vector<DEPEND_VAR> &dv,
                       std::vector<CONS_VAR> &diss, std::vector<CONS_VAR> &rhs,
                       std::vector<PRIM_VAR> &lim, std::vector<PRIM_VAR> &gradx,
                       std::vector<PRIM_VAR> &grady)
    : BaseNumeric(param, geom, cv, dv, diss, rhs, lim, gradx, grady) {}

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

    rrho = 1.0 / cv[i].dens;
    gam1 = dv[i].gamma - 1.0;
    ggm1 = dv[i].gamma / gam1;
    rl = cv[i].dens + lim[i].dens * (gradx[i].dens * rx + grady[i].dens * ry);
    ul = cv[i].xmom * rrho +
         lim[i].velx * (gradx[i].velx * rx + grady[i].velx * ry);
    vl = cv[i].ymom * rrho +
         lim[i].vely * (gradx[i].vely * rx + grady[i].vely * ry);
    pl = dv[i].press +
         lim[i].press * (gradx[i].press * rx + grady[i].press * ry);
    hl = ggm1 * pl / rl + 0.5 * (ul * ul + vl * vl);

    rrho = 1.0 / cv[j].dens;
    gam1 = dv[j].gamma - 1.0;
    ggm1 = dv[j].gamma / gam1;
    rr = cv[j].dens - lim[j].dens * (gradx[j].dens * rx + grady[j].dens * ry);
    ur = cv[j].xmom * rrho -
         lim[j].velx * (gradx[j].velx * rx + grady[j].velx * ry);
    vr = cv[j].ymom * rrho -
         lim[j].vely * (gradx[j].vely * rx + grady[j].vely * ry);
    pr = dv[j].press -
         lim[j].press * (gradx[j].press * rx + grady[j].press * ry);
    hr = ggm1 * pr / rr + 0.5 * (ur * ur + vr * vr);

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
    ry = 0.5 * (geom.coords[i].y - geom.coords[i].y);

    rrho = 1.0 / cv[i].dens;
    gam1 = dv[i].gamma - 1.0;
    ggm1 = dv[i].gamma / gam1;
    rl = cv[i].dens + lim[i].dens * (gradx[i].dens * rx + grady[i].dens * ry);
    ul = cv[i].xmom * rrho +
         lim[i].velx * (gradx[i].velx * rx + grady[i].velx * ry);
    vl = cv[i].ymom * rrho +
         lim[i].vely * (gradx[i].vely * rx + grady[i].vely * ry);
    pl = dv[i].press +
         lim[i].press * (gradx[i].press * rx + grady[i].press * ry);
    hl = ggm1 * pl / rl + 0.5 * (ul * ul + vl * vl);
    qsrl = (ul * geom.sij[ie].x + vl * geom.sij[ie].y) * rl;

    rrho = 1.0 / cv[j].dens;
    gam1 = dv[j].gamma - 1.0;
    ggm1 = dv[j].gamma / gam1;
    rr = cv[j].dens - lim[j].dens * (gradx[j].dens * rx + grady[j].dens * ry);
    ur = cv[j].xmom * rrho -
         lim[j].velx * (gradx[j].velx * rx + grady[j].velx * ry);
    vr = cv[j].ymom * rrho -
         lim[j].vely * (gradx[j].vely * rx + grady[j].vely * ry);
    pr = dv[j].press -
         lim[j].press * (gradx[j].press * rx + grady[j].press * ry);
    hr = ggm1 * pr / rr + 0.5 * (ur * ur + vr * vr);
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
} // namespace solver