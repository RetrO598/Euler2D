#include <solver/FVMSolver.h>

#include <cmath>

namespace solver {
void FVMSolver::ZeroRes() {
  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    if (geom.BoundTypes[ib] >= 500 && geom.BoundTypes[ib] < 600) {
      if (geom.BoundTypes[ib] - 500 < 2) {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.boundaryNode[ibn].node;
          rhs[i].xmom = 0.0;
        }
      } else {
        for (int ibn = ibegn; ibn <= iendn; ++ibn) {
          int i = geom.boundaryNode[ibn].node;
          rhs[i].ymom = 0.0;
        }
      }
    } else if ((geom.BoundTypes[ib] >= 300 && geom.BoundTypes[ib] < 400) &&
               param.equationtype_ == preprocess::equationType::NavierStokes) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.boundaryNode[ib].node;
        rhs[i].xmom = 0.0;
        rhs[i].ymom = 0.0;
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::Irsmoo() {
  for (int i = 0; i < geom.totNodes; ++i) {
    nContr[i] = 0;
    rhsOld[i] = rhs[i];
  }

  for (int itirs = 1; itirs <= param.numOfIterSmooth; ++itirs) {
    for (int i = 0; i < geom.phyNodes; ++i) {
      rhsIter[i].dens = 0.0;
      rhsIter[i].xmom = 0.0;
      rhsIter[i].ymom = 0.0;
      rhsIter[i].ener = 0.0;
    }
    if (itirs == 1) {
      for (int ie = 0; ie < geom.phyEdges; ++ie) {
        int i = geom.edge[ie].nodei;
        int j = geom.edge[ie].nodej;
        nContr[i]++;
        nContr[j]++;

        rhsIter[i].dens += rhs[j].dens;
        rhsIter[i].xmom += rhs[j].xmom;
        rhsIter[i].ymom += rhs[j].ymom;
        rhsIter[i].ener += rhs[j].ener;

        rhsIter[j].dens += rhs[i].dens;
        rhsIter[j].xmom += rhs[i].xmom;
        rhsIter[j].ymom += rhs[i].ymom;
        rhsIter[j].ener += rhs[i].ener;
      }

      PeriodicInt(nContr);
    } else {
      for (int ie = 0; ie < geom.phyEdges; ++ie) {
        int i = geom.edge[ie].nodei;
        int j = geom.edge[ie].nodej;
        rhsIter[i].dens += rhs[j].dens;
        rhsIter[i].xmom += rhs[j].xmom;
        rhsIter[i].ymom += rhs[j].ymom;
        rhsIter[i].ener += rhs[j].ener;

        rhsIter[j].dens += rhs[i].dens;
        rhsIter[j].xmom += rhs[i].xmom;
        rhsIter[j].ymom += rhs[i].ymom;
        rhsIter[j].ener += rhs[i].ener;
      }
    }
    PeriodicCons(rhsIter);

    for (int i = 0; i < geom.phyNodes; ++i) {
      double den = 1.0 / (1.0 + ((double)nContr[i]) * param.imResiSmooth);
      rhs[i].dens =
          (rhsIter[i].dens * param.imResiSmooth + rhsOld[i].dens) * den;
      rhs[i].xmom =
          (rhsIter[i].xmom * param.imResiSmooth + rhsOld[i].xmom) * den;
      rhs[i].ymom =
          (rhsIter[i].ymom * param.imResiSmooth + rhsOld[i].ymom) * den;
      rhs[i].ener =
          (rhsIter[i].ener * param.imResiSmooth + rhsOld[i].ener) * den;
    }
  }
}
} // namespace solver