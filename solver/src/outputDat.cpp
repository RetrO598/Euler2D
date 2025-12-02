#include <solver/FVMSolver.h>

#include <cmath>

#include "element/utils.hpp"
#include "pre/parameter.h"

namespace solver {
void FVMSolver::writeTecplotDat() {
  std::ofstream outFile(param.title + ".dat");
  if (geom.elements[0].getVTKtype() == preprocess::TRIANGLE) {
    outFile << "TITLE = \"Volume Data\"\n";
    outFile << "VARIABLES = \"X\", \"Y\", \"Density\", \"U\", \"V\", \"P\", "
               "\"T\"\n";
    outFile << "ZONE N=" << geom.phyNodes << ", E=" << geom.numElems
            << ", F=FEPOINT, ET=TRIANGLE\n";

    for (int i = 0; i < geom.phyNodes; ++i) {
      outFile << geom.coords[i].x << " " << geom.coords[i].y << " "
              << cv[i].dens << " " << cv[i].xmom / cv[i].dens << " "
              << cv[i].ymom / cv[i].dens << " " << dv[i].press << " "
              << dv[i].temp << "\n";
    }

    for (int i = 0; i < geom.numElems; ++i) {
      outFile << geom.elements[i].getNodes(0) + 1 << " "
              << geom.elements[i].getNodes(1) + 1 << " "
              << geom.elements[i].getNodes(2) + 1 << "\n";
    }
  } else if (geom.elements[0].getVTKtype() == preprocess::RECTANGLE) {
    outFile << "TITLE = \"Volume Data\"\n";
    outFile << "VARIABLES = \"X\", \"Y\", \"Density\", \"U\", \"V\", \"P\", "
               "\"T\"\n";
    outFile << "ZONE N=" << geom.phyNodes << ", E=" << geom.numElems
            << ", F=FEPOINT, ET=QUADRILATERAL\n";

    for (int i = 0; i < geom.phyNodes; ++i) {
      outFile << geom.coords[i].x << " " << geom.coords[i].y << " "
              << cv[i].dens << " " << cv[i].xmom / cv[i].dens << " "
              << cv[i].ymom / cv[i].dens << " " << dv[i].press << " "
              << dv[i].temp << "\n";
    }

    for (int i = 0; i < geom.numElems; ++i) {
      outFile << geom.elements[i].getNodes(0) + 1 << " "
              << geom.elements[i].getNodes(1) + 1 << " "
              << geom.elements[i].getNodes(2) + 1 << " "
              << geom.elements[i].getNodes(3) + 1 << " " << "\n";
    }
  }
  outFile << "\n";
  outFile.close();
}

void FVMSolver::writeLineDat() {
  std::ofstream outFile(param.title + "_line" + ".dat");
  outFile << "TITLE = \"Line Data\"\n";
  outFile
      << "VARIABLES = \"X\", \"Y\", \"Density\", \"U\", \"V\", \"P\", \"T\"\n";
  int nodes = 0;
  int number = 0;
  for (int i = 0; i < geom.numBoundSegs; ++i) {
    auto name = geom.bname[i];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::EulerWall ||
        type == preprocess::BoundaryType::NoSlipWall) {
      nodes += geom.ibound[i].bnodeIndex + 1;
      nodes -= number;
    }
    number = geom.ibound[i].bnodeIndex + 1;
  }
  outFile << "ZONE T=\"Line Data\", I=" << nodes << ", J=1, F=POINT\n";
  int ibeg = 0;
  int iend = 0;
  for (int i = 0; i < geom.numBoundSegs; ++i) {
    iend = geom.ibound[i].bnodeIndex;
    auto name = geom.bname[i];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::EulerWall ||
        type == preprocess::BoundaryType::NoSlipWall) {
      for (int in = ibeg; in <= iend; ++in) {
        outFile << geom.coords[in].x << " " << geom.coords[in].y << " "
                << cv[in].dens << " " << cv[in].xmom / cv[in].dens << " "
                << cv[in].ymom / cv[in].dens << " " << dv[in].press << " "
                << dv[in].temp << "\n";
      }
    }
    ibeg = iend + 1;
  }

  outFile << "\n";
  outFile.close();
}
}  // namespace solver