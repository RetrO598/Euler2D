#include <solver/FVMSolver.h>

#include <cmath>

#include "pre/parameter.h"

namespace solver {
void FVMSolver::Gradients() {
  double fcx[4];
  double fcy[4];

  PRIM_VAR init;
  init.dens = 0.0;
  init.velx = 0.0;
  init.vely = 0.0;
  init.press = 0.0;

  gradx.assign(geom.phyNodes, init);
  grady.assign(geom.phyNodes, init);

  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;

    double rav = 0.5 * (cv[i].dens + cv[j].dens);
    double uav = 0.5 * (cv[i].xmom / cv[i].dens + cv[j].xmom / cv[j].dens);
    double vav = 0.5 * (cv[i].ymom / cv[i].dens + cv[j].ymom / cv[j].dens);
    double pav = 0.5 * (dv[i].press + dv[j].press);

    double sx = geom.sij[ie].x;
    double sy = geom.sij[ie].y;

    fcx[0] = rav * sx;
    fcx[1] = uav * sx;
    fcx[2] = vav * sx;
    fcx[3] = pav * sx;

    fcy[0] = rav * sy;
    fcy[1] = uav * sy;
    fcy[2] = vav * sy;
    fcy[3] = pav * sy;

    gradx[i].dens += fcx[0];
    gradx[i].velx += fcx[1];
    gradx[i].vely += fcx[2];
    gradx[i].press += fcx[3];
    gradx[j].dens -= fcx[0];
    gradx[j].velx -= fcx[1];
    gradx[j].vely -= fcx[2];
    gradx[j].press -= fcx[3];

    grady[i].dens += fcy[0];
    grady[i].velx += fcy[1];
    grady[i].vely += fcy[2];
    grady[i].press += fcy[3];
    grady[j].dens -= fcy[0];
    grady[j].velx -= fcy[1];
    grady[j].vely -= fcy[2];
    grady[j].press -= fcy[3];
  }

  int ibegf = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendf = geom.ibound[ib].bfaceIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    bool flag = true;
    if (type == preprocess::BoundaryType::Symmetric) {
      flag = false;
    }
    if (type == preprocess::BoundaryType::Periodic) {
      flag = false;
    }

    if (flag) {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        int i = geom.boundaryFace[ibf].nodei;
        int j = geom.boundaryFace[ibf].nodej;
        double sx = geom.sbf[ibf].x / 12.0;
        double sy = geom.sbf[ibf].y / 12.0;

        double rav = 5.0 * cv[i].dens + cv[j].dens;
        double uav = 5.0 * cv[i].xmom / cv[i].dens + cv[j].xmom / cv[j].dens;
        double vav = 5.0 * cv[i].ymom / cv[i].dens + cv[j].ymom / cv[j].dens;
        double pav = 5.0 * dv[i].press + dv[j].press;

        gradx[i].dens += rav * sx;
        gradx[i].velx += uav * sx;
        gradx[i].vely += vav * sx;
        gradx[i].press += pav * sx;

        grady[i].dens += rav * sy;
        grady[i].velx += uav * sy;
        grady[i].vely += vav * sy;
        grady[i].press += pav * sy;

        rav = 5.0 * cv[j].dens + cv[i].dens;
        uav = 5.0 * cv[j].xmom / cv[j].dens + cv[i].xmom / cv[i].dens;
        vav = 5.0 * cv[j].ymom / cv[j].dens + cv[i].ymom / cv[i].dens;
        pav = 5.0 * dv[j].press + dv[i].press;

        gradx[j].dens += rav * sx;
        gradx[j].velx += uav * sx;
        gradx[j].vely += vav * sx;
        gradx[j].press += pav * sx;

        grady[j].dens += rav * sy;
        grady[j].velx += uav * sy;
        grady[j].vely += vav * sy;
        grady[j].press += pav * sy;
      }
    }
    ibegf = iendf + 1;
  }

  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Symmetric) {
      if ((geom.BoundTypes[ib] - 500) < 2) {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.boundaryNode[ibn].node;
          gradx[i].dens = 0.0;
          grady[i].velx = 0.0;
          gradx[i].vely = 0.0;
          gradx[i].press = 0.0;
        }
      } else {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.boundaryNode[ibn].node;
          grady[i].dens = 0.0;
          grady[i].velx = 0.0;
          gradx[i].vely = 0.0;
          grady[i].press = 0.0;
        }
      }
    }
    ibegn = iendn + 1;
  }

  PeriodicPrim(gradx);
  PeriodicPrim(grady);

  for (int i = 0; i < geom.phyNodes; ++i) {
    gradx[i].dens /= geom.vol[i];
    gradx[i].velx /= geom.vol[i];
    gradx[i].vely /= geom.vol[i];
    gradx[i].press /= geom.vol[i];

    grady[i].dens /= geom.vol[i];
    grady[i].velx /= geom.vol[i];
    grady[i].vely /= geom.vol[i];
    grady[i].press /= geom.vol[i];
  }
}

void FVMSolver::GradientsVisc() {
  double fcx[5];
  double fcy[5];

  PRIM_VAR init;
  init.dens = 0.0;
  init.velx = 0.0;
  init.vely = 0.0;
  init.press = 0.0;

  gradx.assign(geom.phyNodes, init);
  grady.assign(geom.phyNodes, init);
  gradTx.assign(geom.phyNodes, 0.0);
  gradTy.assign(geom.phyNodes, 0.0);

  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;

    double rav = 0.5 * (cv[i].dens + cv[j].dens);
    double uav = 0.5 * (cv[i].xmom / cv[i].dens + cv[j].xmom / cv[j].dens);
    double vav = 0.5 * (cv[i].ymom / cv[i].dens + cv[j].ymom / cv[j].dens);
    double pav = 0.5 * (dv[i].press + dv[j].press);
    double tav = 0.5 * (dv[i].temp + dv[j].temp);

    double sx = geom.sij[ie].x;
    double sy = geom.sij[ie].y;

    fcx[0] = rav * sx;
    fcx[1] = uav * sx;
    fcx[2] = vav * sx;
    fcx[3] = pav * sx;
    fcx[4] = tav * sx;

    fcy[0] = rav * sy;
    fcy[1] = uav * sy;
    fcy[2] = vav * sy;
    fcy[3] = pav * sy;
    fcy[4] = tav * sy;

    gradx[i].dens += fcx[0];
    gradx[i].velx += fcx[1];
    gradx[i].vely += fcx[2];
    gradx[i].press += fcx[3];
    gradTx[i] += fcx[4];
    gradx[j].dens -= fcx[0];
    gradx[j].velx -= fcx[1];
    gradx[j].vely -= fcx[2];
    gradx[j].press -= fcx[3];
    gradTx[j] -= fcx[4];

    grady[i].dens += fcy[0];
    grady[i].velx += fcy[1];
    grady[i].vely += fcy[2];
    grady[i].press += fcy[3];
    gradTy[i] += fcy[4];
    grady[j].dens -= fcy[0];
    grady[j].velx -= fcy[1];
    grady[j].vely -= fcy[2];
    grady[j].press -= fcy[3];
    gradTy[j] -= fcy[4];
  }

  int ibegf = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendf = geom.ibound[ib].bfaceIndex;
    bool flag = true;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Symmetric) {
      flag = false;
    }
    if (type == preprocess::BoundaryType::Periodic) {
      flag = false;
    }

    if (flag) {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        int i = geom.boundaryFace[ibf].nodei;
        int j = geom.boundaryFace[ibf].nodej;
        double sx = geom.sbf[ibf].x / 12.0;
        double sy = geom.sbf[ibf].y / 12.0;

        double rav = 5.0 * cv[i].dens + cv[j].dens;
        double uav = 5.0 * cv[i].xmom / cv[i].dens + cv[j].xmom / cv[j].dens;
        double vav = 5.0 * cv[i].ymom / cv[i].dens + cv[j].ymom / cv[j].dens;
        double pav = 5.0 * dv[i].press + dv[j].press;
        double tav = 5.0 * dv[i].temp + dv[j].temp;

        gradx[i].dens += rav * sx;
        gradx[i].velx += uav * sx;
        gradx[i].vely += vav * sx;
        gradx[i].press += pav * sx;
        gradTx[i] += tav * sx;

        grady[i].dens += rav * sy;
        grady[i].velx += uav * sy;
        grady[i].vely += vav * sy;
        grady[i].press += pav * sy;
        gradTy[i] += tav * sy;

        rav = 5.0 * cv[j].dens + cv[i].dens;
        uav = 5.0 * cv[j].xmom / cv[j].dens + cv[i].xmom / cv[i].dens;
        vav = 5.0 * cv[j].ymom / cv[j].dens + cv[i].ymom / cv[i].dens;
        pav = 5.0 * dv[j].press + dv[i].press;
        tav = 5.0 * dv[j].temp + dv[i].temp;

        gradx[j].dens += rav * sx;
        gradx[j].velx += uav * sx;
        gradx[j].vely += vav * sx;
        gradx[j].press += pav * sx;
        gradTx[j] += tav * sx;

        grady[j].dens += rav * sy;
        grady[j].velx += uav * sy;
        grady[j].vely += vav * sy;
        grady[j].press += pav * sy;
        gradTy[j] += tav * sy;
      }
    }
    ibegf = iendf + 1;
  }

  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Symmetric) {
      if ((geom.BoundTypes[ib] - 500) < 2) {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.boundaryNode[ibn].node;
          gradx[i].dens = 0.0;
          grady[i].velx = 0.0;
          gradx[i].vely = 0.0;
          gradx[i].press = 0.0;
          gradTx[i] = 0.0;
        }
      } else {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.boundaryNode[ibn].node;
          grady[i].dens = 0.0;
          grady[i].velx = 0.0;
          gradx[i].vely = 0.0;
          grady[i].press = 0.0;
          gradTy[i] = 0.0;
        }
      }
    }
    ibegn = iendn + 1;
  }

  PeriodicPrim(gradx);
  PeriodicPrim(grady);
  PeriodicVisc(gradTx, gradTy);

  for (int i = 0; i < geom.phyNodes; ++i) {
    gradx[i].dens /= geom.vol[i];
    gradx[i].velx /= geom.vol[i];
    gradx[i].vely /= geom.vol[i];
    gradx[i].press /= geom.vol[i];
    gradTx[i] /= geom.vol[i];

    grady[i].dens /= geom.vol[i];
    grady[i].velx /= geom.vol[i];
    grady[i].vely /= geom.vol[i];
    grady[i].press /= geom.vol[i];
    gradTy[i] /= geom.vol[i];
  }
}
}  // namespace solver