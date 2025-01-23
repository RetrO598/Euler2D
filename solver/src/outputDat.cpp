#include <solver/FVMSolver.h>

#include <cmath>

namespace solver {
void FVMSolver::writeTecplotDat() {
  std::ofstream outFile(param.title + ".dat");
  outFile << "TITLE = \"Unstructured Triangle Mesh\"\n";
  outFile
      << "VARIABLES = \"X\", \"Y\", \"Density\", \"U\", \"V\", \"P\", \"T\"\n";
  outFile << "ZONE N=" << geom.phyNodes << ", E=" << geom.numTria
          << ", F=FEPOINT, ET=TRIANGLE\n";

  for (int i = 0; i < geom.phyNodes; ++i) {
    outFile << geom.coords[i].x << " " << geom.coords[i].y << " " << cv[i].dens
            << " " << cv[i].xmom / cv[i].dens << " " << cv[i].ymom / cv[i].dens
            << " " << dv[i].press << " " << dv[i].temp << "\n";
  }

  for (int i = 0; i < geom.numTria; ++i) {
    outFile << geom.tria[i].node[0] + 1 << " " << geom.tria[i].node[1] + 1
            << " " << geom.tria[i].node[2] + 1 << "\n";
  }
  outFile << "\n";
  outFile.close();
}
}  // namespace solver