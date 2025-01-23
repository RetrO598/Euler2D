#include <solver/FVMSolver.h>

#include <cmath>

namespace solver {
void FVMSolver::BoundaryConditions() {
  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int itype = geom.BoundTypes[ib];
    int iendn = geom.ibound[ib].bnodeIndex;

    if (itype >= 100 && itype < 200) {
      BoundInflow(ibegn, iendn);
    } else if (itype >= 200 && itype < 300) {
      BoundOutflow(ibegn, iendn);
    } else if (itype >= 600 && itype < 700) {
      BoundFarfield(ibegn, iendn);
    }
    ibegn = iendn + 1;
  }

  ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int itype = geom.BoundTypes[ib];
    int iendn = geom.ibound[ib].bnodeIndex;
    if (param.equationtype_ == preprocess::equationType::NavierStokes &&
        (itype >= 300 && itype < 400)) {
      BoundWallVisc(ibegn, iendn);
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::BoundInflow(int beg, int end) {
  for (int ib = beg; ib <= end; ++ib) {
    int ibn = geom.boundNode[ib].node;
    int idn = geom.boundNode[ib].dummy;
    int ie = geom.boundNode[ib].indexEdge;

    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double sxn = geom.sij[ie].x / ds;
    double syn = geom.sij[ie].y / ds;

    double gam1 = dv[ibn].gamma - 1.0;
    double ggm1 = dv[ibn].gamma / gam1;
    double rgas = gam1 * dv[ibn].Cp / dv[ibn].gamma;
    double rrho = 1.0 / cv[ibn].dens;
    double u = cv[ibn].xmom * rrho;
    double v = cv[ibn].ymom * rrho;
    double uabs = std::sqrt(u * u + v * v);
    double unorm = u * sxn + v * syn;
    double cosa = 0.0;
    if (uabs < 1e-16) {
      cosa = 1.0;
    } else {
      cosa = -unorm / uabs;
    }

    double c02 = dv[ibn].cs * dv[ibn].cs + 0.5 * gam1 * uabs * uabs;
    double rinv = unorm - 2.0 * dv[ibn].cs / gam1;
    double dis =
        (gam1 * cosa * cosa + 2.0) * c02 / (gam1 * rinv * rinv) - 0.5 * gam1;
    dis = std::max(dis, 1e-16);
    double cb = -rinv * (gam1 / (gam1 * cosa * cosa + 2.0)) *
                (1.0 + cosa * std::sqrt(dis));
    double cc02 = std::min((cb * cb) / c02, 1.0);
    double tb = param.TtInlet * cc02;
    double pb = param.PtInlet * std::pow((tb / param.TtInlet), ggm1);
    double rhob = pb / (rgas * tb);
    double uabsb = (2.0 * dv[ibn].Cp) * (param.TtInlet - tb);
    uabsb = std::sqrt(uabsb);
    double ub = uabsb * std::cos(param.flowAngIn);
    double vb = uabsb * std::sin(param.flowAngIn);

    cv[idn].dens = rhob;
    cv[idn].xmom = rhob * ub;
    cv[idn].ymom = rhob * vb;
    cv[idn].ener = pb / gam1 + 0.5 * rhob * (ub * ub + vb * vb);

    ConvToDepend(idn);
  }
}

void FVMSolver::BoundFarfield(int beg, int end) {
  PRIM_VAR farfield;
  if (param.farCorrect) {
    double bet = std::sqrt(1.0 - param.MaInfinity * param.MaInfinity);
    double cir = 0.25 * param.chord * cl * param.velInfinity / M_PI;
    for (int ib = beg; ib <= end; ++ib) {
      int ibn = geom.boundNode[ib].node;
      int idn = geom.boundNode[ib].dummy;
      int ie = geom.boundNode[ib].indexEdge;

      double gam1 = dv[ibn].gamma - 1.0;
      double ggm1 = dv[ibn].gamma / gam1;
      double gmr = 1.0 / dv[ibn].gamma;
      double gmg = gam1 / dv[ibn].gamma;
      double xa = geom.coords[ibn].x - param.xRefPoint;
      double ya = geom.coords[ibn].y - param.yRefPoint;
      double dist = std::sqrt(xa * xa + ya * ya);
      double angle = std::atan2(ya, xa);
      double sn = std::sin(angle - param.Aoa);
      double dn = 1.0 - param.MaInfinity * param.MaInfinity * sn * sn;
      double vc = cir * bet / (dn * dist);

      farfield.velx = param.uInfinity + vc * std::sin(angle);
      farfield.vely = param.vInfinity - vc * std::cos(angle);
      double qv2 =
          farfield.velx * farfield.velx + farfield.vely * farfield.vely;
      farfield.press =
          std::pow(std::pow(param.PsInfinity, gmg) +
                       gmg * param.RhoInfinity *
                           (param.velInfinity * param.velInfinity - qv2) /
                           (2.0 * std::pow(param.PsInfinity, gmr)),
                   ggm1);
      farfield.dens =
          param.RhoInfinity * std::pow(farfield.press / param.PsInfinity, gmr);

      double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                            geom.sij[ie].y * geom.sij[ie].y);
      double sxn = geom.sij[ie].x / ds;
      double syn = geom.sij[ie].y / ds;

      double rhoe = cv[ibn].dens;
      double ue = cv[ibn].xmom / rhoe;
      double ve = cv[ibn].ymom / rhoe;
      double pe = dv[ibn].press;
      double qn = sxn * ue + syn * ve;
      double crho0 = dv[ibn].cs * rhoe;
      double rhoa;
      double ua;
      double va;
      double pa;
      double sgn;
      double pb;

      if (param.MaInfinity < 1.0) {
        if (qn < 0.0) {
          rhoa = farfield.dens;
          ua = farfield.velx;
          va = farfield.vely;
          pa = farfield.press;
          sgn = -1.0;
          pb = 0.5 * (pa + pe - crho0 * (sxn * (ua - ue) + syn * (va - ve)));
        } else {
          rhoa = rhoe;
          ua = ue;
          va = ve;
          pa = pe;
          sgn = 1.0;
          pb = farfield.press;
        }
        double dprhoc = sgn * (pa - pb) / crho0;
        cv[idn].dens = rhoa + (pb - pa) / (dv[ibn].cs * dv[ibn].cs);
        cv[idn].xmom = cv[idn].dens * (ua + sxn * dprhoc);
        cv[idn].ymom = cv[idn].dens * (va + syn * dprhoc);
        cv[idn].ener =
            pb / gam1 +
            0.5 * (cv[idn].xmom * cv[idn].xmom + cv[idn].ymom * cv[idn].ymom) /
                cv[idn].dens;
      } else {
        if (qn < 0.0) {
          cv[idn].dens = param.RhoInfinity;
          cv[idn].xmom = param.RhoInfinity * param.uInfinity;
          cv[idn].ymom = param.RhoInfinity * param.vInfinity;
          cv[idn].ener = param.PsInfinity / gam1 + 0.5 * param.RhoInfinity *
                                                       param.velInfinity *
                                                       param.velInfinity;
        } else {
          cv[idn].dens = cv[ibn].dens;
          cv[idn].xmom = cv[ibn].xmom;
          cv[idn].ymom = cv[ibn].ymom;
          cv[idn].ener = cv[ibn].ener;
        }
      }
      ConvToDepend(idn);
    }
  } else {
    for (int ib = beg; ib <= end; ++ib) {
      int ibn = geom.boundNode[ib].node;
      int idn = geom.boundNode[ib].dummy;
      int ie = geom.boundNode[ib].indexEdge;

      double gam1 = dv[ibn].gamma - 1.0;
      farfield.dens = param.RhoInfinity;
      farfield.velx = param.uInfinity;
      farfield.vely = param.vInfinity;
      farfield.press = param.PsInfinity;

      double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                            geom.sij[ie].y * geom.sij[ie].y);
      double sxn = geom.sij[ie].x / ds;
      double syn = geom.sij[ie].y / ds;

      double rhoe = cv[ibn].dens;
      double ue = cv[ibn].xmom / rhoe;
      double ve = cv[ibn].ymom / rhoe;
      double pe = dv[ibn].press;
      double qn = sxn * ue + syn * ve;
      double crho0 = dv[ibn].cs * rhoe;
      double rhoa;
      double ua;
      double va;
      double pa;
      double sgn;
      double pb;

      if (param.MaInfinity < 1.0) {
        if (qn < 0.0) {
          rhoa = farfield.dens;
          ua = farfield.velx;
          va = farfield.vely;
          pa = farfield.press;
          sgn = -1.0;
          pb = 0.5 * (pa + pe - crho0 * (sxn * (ua - ue) + syn * (va - ve)));
        } else {
          rhoa = rhoe;
          ua = ue;
          va = ve;
          pa = pe;
          sgn = 1.0;
          pb = farfield.press;
        }
        double dprhoc = sgn * (pa - pb) / crho0;
        cv[idn].dens = rhoa + (pb - pa) / (dv[ibn].cs * dv[ibn].cs);
        cv[idn].xmom = cv[idn].dens * (ua + sxn * dprhoc);
        cv[idn].ymom = cv[idn].dens * (va + syn * dprhoc);
        cv[idn].ener =
            pb / gam1 +
            0.5 * (cv[idn].xmom * cv[idn].xmom + cv[idn].ymom * cv[idn].ymom) /
                cv[idn].dens;
      } else {
        if (qn < 0.0) {
          cv[idn].dens = param.RhoInfinity;
          cv[idn].xmom = param.RhoInfinity * param.uInfinity;
          cv[idn].ymom = param.RhoInfinity * param.vInfinity;
          cv[idn].ener = param.PsInfinity / gam1 + 0.5 * param.RhoInfinity *
                                                       param.velInfinity *
                                                       param.velInfinity;
        } else {
          cv[idn].dens = cv[ibn].dens;
          cv[idn].xmom = cv[ibn].xmom;
          cv[idn].ymom = cv[ibn].ymom;
          cv[idn].ener = cv[ibn].ener;
        }
      }
      ConvToDepend(idn);
    }
  }
}

void FVMSolver::BoundOutflow(int beg, int end) {
  for (int ib = beg; ib <= end; ++ib) {
    int ibn = geom.boundNode[ib].node;
    int idn = geom.boundNode[ib].dummy;
    int ie = geom.boundNode[ib].indexEdge;

    double ds = std::sqrt(geom.sij[ie].x * geom.sij[ie].x +
                          geom.sij[ie].y * geom.sij[ie].y);
    double sxn = geom.sij[ie].x / ds;
    double syn = geom.sij[ie].y / ds;

    double gam1 = dv[ibn].gamma - 1.0;
    double rrho = 1.0 / cv[ibn].dens;
    double u = cv[ibn].xmom * rrho;
    double v = cv[ibn].ymom * rrho;
    double q = std::sqrt(u * u + v * v);
    double mach = q / dv[ibn].cs;

    if (mach < 1.0) {
      double rrhoc = rrho / dv[ibn].cs;
      double deltp = dv[ibn].press - param.PsOutlet;
      double rhob = cv[ibn].dens - deltp / (dv[ibn].cs * dv[ibn].cs);
      double ub = u + sxn * deltp * rrhoc;
      double vb = v + syn * deltp * rrhoc;
      double vnd = ub * sxn + vb * syn;

      if (vnd < 0.0) {
        ub = (((u) < (0)) ? (-1.0) : (1.0)) *
             std::max(std::abs(ub), std::abs(u));
        vb = (((v) < (0)) ? (-1.0) : (1.0)) *
             std::max(std::abs(vb), std::abs(v));
      }
      cv[idn].dens = rhob;
      cv[idn].xmom = rhob * ub;
      cv[idn].ymom = rhob * vb;
      cv[idn].ener = param.PsOutlet / gam1 + 0.5 * rhob * (ub * ub + vb * vb);
    } else {
      cv[idn].dens = cv[ibn].dens;
      cv[idn].xmom = cv[ibn].xmom;
      cv[idn].ymom = cv[ibn].ymom;
      cv[idn].ener = cv[ibn].ener;
    }
    ConvToDepend(idn);
  }
}

void FVMSolver::BoundWallVisc(int beg, int end) {
  for (int ib = beg; ib <= end; ++ib) {
    int ibn = geom.boundNode[ib].node;

    cv[ibn].xmom = 0.0;
    cv[ibn].ymom = 0.0;

    ConvToDepend(ibn);
  }
}

void FVMSolver::PeriodicPrim(std::vector<PRIM_VAR> &var) {
  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    if (geom.BoundTypes[ib] >= 700 && geom.BoundTypes[ib] < 800) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.boundNode[ibn].node;
        int j = geom.boundNode[ibn].dummy;

        var[i].dens += var[j].dens;
        var[i].velx += var[j].velx;
        var[i].vely += var[j].vely;
        var[i].press += var[j].press;
        var[j].dens = var[i].dens;
        var[j].velx = var[i].velx;
        var[j].vely = var[i].vely;
        var[j].press = var[i].press;
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::PeriodicCons(std::vector<CONS_VAR> &var) {
  int ibegn = 0.0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    if (geom.BoundTypes[ib] >= 700 && geom.BoundTypes[ib] < 800) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.boundNode[ibn].node;
        int j = geom.boundNode[ibn].dummy;
        var[i].dens += var[j].dens;
        var[i].xmom += var[j].xmom;
        var[i].ymom += var[j].ymom;
        var[i].ener += var[j].ener;
        var[j].dens = var[i].dens;
        var[j].xmom = var[i].xmom;
        var[j].ymom = var[i].ymom;
        var[j].ener = var[i].ener;
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::PeriodicInt(std::vector<int> &var) {
  int i, j, ib, ibn, ibegn, iendn;

  ibegn = 0;

  for (ib = 0; ib < geom.numBoundSegs; ib++) {
    iendn = geom.ibound[ib].bnodeIndex;
    if (geom.BoundTypes[ib] >= 700 && geom.BoundTypes[ib] < 800) {
      for (ibn = ibegn; ibn <= iendn; ibn++) {
        i = geom.boundNode[ibn].node;
        j = geom.boundNode[ibn].dummy;
        var[i] += var[j];
        var[j] = var[i];
      }
    }
    ibegn = iendn + 1;
  }
}

}  // namespace solver