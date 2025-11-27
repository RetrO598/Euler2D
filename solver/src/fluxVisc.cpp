
#include <solver/FVMSolver.h>

#include <cmath>

namespace solver {
void FVMSolver::fluxVisc(double beta) {
  double fv[3];

  for (int ie = 0; ie < geom.phyEdges; ++ie) {
    int i = geom.edge[ie].nodei;
    int j = geom.edge[ie].nodej;

    double ui = cv[i].xmom / cv[i].dens;
    double uj = cv[j].xmom / cv[j].dens;
    double vi = cv[i].ymom / cv[i].dens;
    double vj = cv[j].ymom / cv[j].dens;
    double uav = 0.5 * (ui + uj);
    double vav = 0.5 * (vi + vj);
    double mav = 0.5 * (dvlam[i].mu + dvlam[j].mu);
    double kav = 0.5 * (dvlam[i].lambda + dvlam[j].lambda);

    double txn = geom.coords[j].x - geom.coords[i].x;
    double tyn = geom.coords[j].y - geom.coords[i].y;
    double ds2 = txn * txn + tyn * tyn;
    double rds = 1.0 / std::sqrt(ds2);
    txn = txn * rds;
    tyn = tyn * rds;

    double duxa = 0.5 * (gradx[i].velx + gradx[j].velx);
    double duya = 0.5 * (grady[i].velx + grady[j].velx);
    double dvxa = 0.5 * (gradx[i].vely + gradx[j].vely);
    double dvya = 0.5 * (grady[i].vely + grady[j].vely);
    double dtxa = 0.5 * (gradTx[i] + gradTx[j]);
    double dtya = 0.5 * (gradTy[i] + gradTy[j]);

    double duds = rds * (uj - ui);
    double dvds = rds * (vj - vi);
    double dtds = rds * (dv[j].temp - dv[i].temp);

    double dudt = duxa * txn + duya * tyn - duds;
    double dvdt = dvxa * txn + dvya * tyn - dvds;
    double dtdt = dtxa * txn + dtya * tyn - dtds;

    double duxf = duxa - dudt * txn;
    double duyf = duya - dudt * tyn;
    double dvxf = dvxa - dvdt * txn;
    double dvyf = dvya - dvdt * tyn;
    double dtxf = dtxa - dtdt * txn;
    double dtyf = dtya - dtdt * tyn;

    double tauxx = 2.0 / 3.0 * mav * (2.0 * duxf - dvyf);
    double tauyy = 2.0 / 3.0 * mav * (2.0 * dvyf - duxf);
    double tauxy = mav * (duyf + dvxf);
    double phix = uav * tauxx + vav * tauxy + kav * dtxf;
    double phiy = uav * tauxy + vav * tauyy + kav * dtyf;

    fv[0] = geom.sij[ie].x * tauxx + geom.sij[ie].y * tauxy;
    fv[1] = geom.sij[ie].x * tauxy + geom.sij[ie].y * tauyy;
    fv[2] = geom.sij[ie].x * phix + geom.sij[ie].y * phiy;

    diss[i].xmom += fv[0];
    diss[i].ymom += fv[1];
    diss[i].ener += fv[2];

    diss[j].xmom -= fv[0];
    diss[j].ymom -= fv[1];
    diss[j].ener -= fv[2];
  }
}
}  // namespace solver