#include <solver/FVMSolver.h>

#include <cmath>
#include <stdexcept>
#include <vector>

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
  outFile << "VARIABLES = \"X\", \"Y\", \"Density\", \"U\", \"V\", \"P\", "
             "\"T\", \"Ma\"\n";
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
  std::vector<int> indexVector;
  indexVector.reserve(nodes);

  int ibegn = 0;
  int ibegf = 0;
  for (int i = 0; i < geom.numBoundSegs; ++i) {
    auto name = geom.bname[i];
    auto type = param.boundaryMap.find(name)->second;
    int iendn = geom.ibound[i].bnodeIndex;
    int iendf = geom.ibound[i].bfaceIndex;
    if (type == preprocess::BoundaryType::EulerWall ||
        type == preprocess::BoundaryType::NoSlipWall) {
      std::unordered_map<int, std::vector<int>> adj;
      for (int iface = ibegf; iface <= iendf; ++iface) {
        adj[geom.boundaryFace[iface].nodei].push_back(
            geom.boundaryFace[iface].nodej);
        adj[geom.boundaryFace[iface].nodej].push_back(
            geom.boundaryFace[iface].nodei);
      }
      if (adj.empty()) {
        continue;
      }

      int start_node = adj.begin()->first;
      for (auto const& [node, neighs] : adj) {
        if (neighs.size() == 1) {
          start_node = node;
          std::cout << start_node << "\n";
          break;
        }
      }

      int prev = -1;
      int cur = start_node;
      while (cur != -1) {
        indexVector.push_back(cur);
        int next = -1;
        for (int v : adj[cur]) {
          if (v != prev) {
            next = v;
            break;
          }
        }
        prev = cur;
        cur = next;
        if (cur == start_node) break;
      }
    }
    ibegn = iendn + 1;
    ibegf = iendf + 1;
  }

  if (indexVector.size() != nodes) {
    throw std::runtime_error("Number of wall boundary vertex not match.\n");
  }

  outFile << "ZONE T=\"Line Data\", I=" << nodes << ", J=1, F=POINT\n";

  for (auto& i : indexVector) {
    outFile << geom.coords[i].x << " " << geom.coords[i].y << " " << cv[i].dens
            << " " << cv[i].xmom / cv[i].dens << " " << cv[i].ymom / cv[i].dens
            << " " << dv[i].press << " " << dv[i].temp << " "
            << std::sqrt(std::pow(cv[i].xmom / cv[i].dens, 2) +
                         std::pow(cv[i].ymom / cv[i].dens, 2)) /
                   dv[i].cs
            << "\n";
  }

  outFile << "\n";
  outFile.close();
}
}  // namespace solver