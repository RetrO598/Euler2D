#include <pre/geometry.h>
#include <pre/macro.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "common/config.hpp"
#include "pre/parameter.h"

namespace preprocess {
Geometry::Geometry(const parameter& param, const std::string& commentChar)
    : boundaryMap(std::move(param.boundaryMap)),
      periodicInfo(param.periodics),
      periodicMaster(param.periodicMaster) {
  phyNodes = 0;
  phyEdges = 0;
  numElems = 0;
  numBoundSegs = 0;
  numBoundFaces = 0;
  numBoundNodes = 0;
}

void Geometry::printInfo() {
  std::cout << "No. of interior nodes: " << phyNodes << "\n";
  std::cout << "No. of grid cells :" << numElems << "\n";
  std::cout << "No. of interior edges: " << phyEdges << "\n";
  std::cout << "Total number of edges: " << phyEdges << "\n";
  std::cout << "No. of boundary faces: " << numBoundFaces << "\n";
  std::cout << "No. of boundary nodes: " << numBoundNodes << "\n";
}

void Geometry::CheckMetrics() {
  int i, j, ib, ibf, ibn, ie, ibegf, iendf, ibegn, iendn;
  double volmin, volmax, s, smax;
  std::vector<Node> fvecSum;

  volmin = +1.0e+32;
  volmax = -1.0e+32;

  for (i = 0; i < phyNodes; ++i) {
    volmin = std::min(volmin, vol[i]);
    volmax = std::max(volmax, vol[i]);
  }

  fvecSum.resize(phyNodes);
  Node node{0.0, 0.0};
  std::fill(fvecSum.begin(), fvecSum.end(), node);

  for (ie = 0; ie < phyEdges; ++ie) {
    i = edge[ie].nodei;
    j = edge[ie].nodej;
    fvecSum[i].x += sij[ie].x;
    fvecSum[i].y += sij[ie].y;
    fvecSum[j].x -= sij[ie].x;
    fvecSum[j].y -= sij[ie].y;
  }

  ibegf = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendf = ibound[ib].bfaceIndex;
    auto name = bname[ib];
    std::cout << name << "\n";
    auto type = boundaryMap.find(name)->second;
    if (type != BoundaryType::Periodic) {
      for (ibf = ibegf; ibf <= iendf; ++ibf) {
        i = boundaryFace[ibf].nodei;
        j = boundaryFace[ibf].nodej;
        fvecSum[i].x += 0.5 * sbf[ibf].x;
        fvecSum[i].y += 0.5 * sbf[ibf].y;
        fvecSum[j].x += 0.5 * sbf[ibf].x;
        fvecSum[j].y += 0.5 * sbf[ibf].y;
      }
    }

    ibegf = iendf + 1;
  }

  ibegn = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendn = ibound[ib].bnodeIndex;
    auto name = bname[ib];
    auto type = boundaryMap.find(name)->second;
    if (type == BoundaryType::Periodic &&
        periodicMaster.find(name) != periodicMaster.end()) {
      for (ibn = ibegn; ibn <= iendn; ++ibn) {
        i = vertexList[ibn].nodeIdx;
        j = vertexList[ibn].periodicPair;
        fvecSum[i].x += fvecSum[j].x;
        fvecSum[i].y += fvecSum[j].y;
        fvecSum[j].x = fvecSum[i].x;
        fvecSum[j].y = fvecSum[i].y;
      }
    }
    ibegn = iendn + 1;
  }

  smax = -1.0e+32;
  for (i = 0; i < phyNodes; ++i) {
    s = std::sqrt(fvecSum[i].x * fvecSum[i].x + fvecSum[i].y * fvecSum[i].y);
    smax = std::max(smax, s);
  }

  std::cout << " max. sum(S) = " << smax << "\n";
  std::cout << " min. volume = " << volmin << "\n";
  std::cout << " max. volume = " << volmax << "\n";
}

void Geometry::FaceVectorsSymm() {
  int i, j, ib, ibf, ibn, ie, ibegf, ibegn, iendf, iendn;
  std::vector<int> marker;
  double sx, sy;

  marker.resize(phyNodes);

  std::fill(marker.begin(), marker.end(), -1);

  ibegf = 0;
  ibegn = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendf = ibound[ib].bfaceIndex;
    iendn = ibound[ib].bnodeIndex;
    auto name = bname[ib];
    auto type = boundaryMap.find(name)->second;
    if (type == BoundaryType::Symmetric) {
      int symmetryIndex;
      sx = 0.0;
      sy = 0.0;
      for (ibf = ibegf; ibf <= iendf; ++ibf) {
        sx += sbf[ibf].x;
        sy += sbf[ibf].y;
      }
      if (std::abs(sx) > std::abs(sy)) {
        symmetryIndex = 501;
      } else {
        symmetryIndex = 502;
      }

      for (ibn = ibegn; ibn <= iendn; ++ibn) {
        marker[vertexList[ibn].nodeIdx] = symmetryIndex - 500;
      }
    }
    ibegf = iendf + 1;
    ibegn = iendn + 1;
  }

  for (ie = 0; ie < phyEdges; ++ie) {
    i = edge[ie].nodei;
    j = edge[ie].nodej;
    if (marker[i] != -1 && marker[j] != -1) {
      if (marker[i] < 2)  // x=const. plane
        sij[ie].x = 0.0;
      else  // y=const. plane
        sij[ie].y = 0.0;
    }
  }
}

void Geometry::volumeProjections() {
  int ibegn, iendn, ibegf, iendf;
  int i, j, ib, ibf, ibn, ie;
  double sx, sy;

  for (i = 0; i < phyNodes; ++i) {
    sproj[i].x = 0.0;
    sproj[i].y = 0.0;
  }

  for (ie = 0; ie < phyEdges; ++ie) {
    i = edge[ie].nodei;
    j = edge[ie].nodej;
    sx = 0.5 * std::fabs(sij[ie].x);
    sy = 0.5 * std::fabs(sij[ie].y);
    sproj[i].x += sx;
    sproj[i].y += sy;
    sproj[j].x += sx;
    sproj[j].y += sy;
  }

  ibegf = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendf = ibound[ib].bfaceIndex;
    auto name = bname[ib];
    auto type = boundaryMap.find(name)->second;
    if (type != BoundaryType::Periodic) {
      for (ibf = ibegf; ibf <= iendf; ++ibf) {
        i = boundaryFace[ibf].nodei;
        j = boundaryFace[ibf].nodej;
        sx = 0.25 * std::fabs(sbf[ibf].x);
        sy = 0.25 * std::fabs(sbf[ibf].y);
        sproj[i].x += sx;
        sproj[i].y += sy;
        sproj[j].x += sx;
        sproj[j].y += sy;
      }
    }
    ibegf = iendf + 1;
  }

  ibegn = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendn = ibound[ib].bnodeIndex;
    auto name = bname[ib];
    auto type = boundaryMap.find(name)->second;
    if (type == BoundaryType::Periodic &&
        periodicMaster.find(name) != periodicMaster.end()) {
      for (ibn = ibegn; ibn <= iendn; ibn++) {
        i = vertexList[ibn].nodeIdx;
        j = vertexList[ibn].periodicPair;
        sproj[i].x += sproj[j].x;
        sproj[i].y += sproj[j].y;
        sproj[j].x = sproj[i].x;
        sproj[j].y = sproj[i].y;
      }
    }
    ibegn = iendn + 1;
  }
}

void Geometry::ComputeWallDistance() {
  Index nWallVertex = 0;
  Index nodeNum = 0;
  for (std::size_t i = 0; i < numBoundSegs; ++i) {
    auto name = bname[i];
    auto type = boundaryMap.find(name)->second;
    if (type == BoundaryType::NoSlipWall) {
      nWallVertex += ibound[i].bnodeIndex + 1;
      nWallVertex -= nodeNum;
    }
    nodeNum = ibound[i].bnodeIndex + 1;
  }

  std::vector<std::array<double, 2>> coords;
  coords.reserve(nWallVertex);

  nWallVertex = 0;
  nodeNum = 0;
  Index ibegn = 0;
  Index iendn = 0;
  for (std::size_t i = 0; i < numBoundSegs; ++i) {
    auto name = bname[i];
    auto type = boundaryMap.find(name)->second;
    iendn = ibound[i].bnodeIndex;
    if (type == BoundaryType::NoSlipWall) {
      nWallVertex += ibound[i].bnodeIndex + 1;
      nWallVertex -= nodeNum;
      for (Index inode = ibegn; inode <= iendn; ++inode) {
        Index nodeIndex = vertexList[inode].nodeIdx;
        coords.push_back({pointList[nodeIndex].getCoord(0),
                          pointList[nodeIndex].getCoord(1)});
      }
    }
    nodeNum = ibound[i].bnodeIndex + 1;
    ibegn = iendn + 1;
  }

  if (coords.size() != nWallVertex) {
    throw std::runtime_error("number of wall vertex not matching\n");
  }

  if (nWallVertex != 0) {
    for (auto& ipoint : pointList) {
      double distance = std::numeric_limits<double>::max();
      std::array<double, 2> coord{ipoint.getCoord(0), ipoint.getCoord(1)};
      double dist2 = 0.0;
      for (auto& coordBound : coords) {
        double dx = coord[0] - coordBound[0];
        double dy = coord[1] - coordBound[1];
        dist2 = dx * dx + dy * dy;
        distance = std::min(dist2, distance);
      }
      ipoint.setWallDistance(std::sqrt(distance));
    }
  } else {
    for (auto& ipoint : pointList) {
      ipoint.setWallDistance(0.0);
    }
  }
}

void Geometry::MatchPeriodic() {
  int ibegn = 0;
  for (std::size_t i = 0; i < numBoundSegs; ++i) {
    auto name = bname[i];
    auto type = boundaryMap.find(name)->second;
    int iendn = ibound[i].bnodeIndex;
    if (type == preprocess::BoundaryType::Periodic &&
        periodicMaster.find(name) != periodicMaster.end()) {
      std::string slave;
      std::array<double, 2> rotate_center;
      double rotation;
      std::array<double, 2> translation;
      for (auto& info : periodicInfo) {
        if (info.master_name == name) {
          slave = info.slave_name;
          rotate_center = info.center;
          rotation = info.rotate;
          translation = info.translate;
        }
      }

      std::size_t slave_index;
      for (std::size_t iname = 0; iname < bname.size(); ++iname) {
        if (bname[iname] == slave) {
          slave_index = iname;
        }
      }

      std::size_t slave_begn = ibound[slave_index - 1].bnodeIndex + 1;
      std::size_t slave_endn = ibound[slave_index].bnodeIndex;

      translation[0] += rotate_center[0];
      translation[1] += rotate_center[1];

      double Psi = rotation;
      double sinPsi = std::sin(Psi);
      double cosPsi = std::cos(Psi);

      for (std::size_t inode = ibegn; inode <= iendn; ++inode) {
        std::size_t ibn = vertexList[inode].nodeIdx;

        auto coordx = coords[ibn].x;
        auto coordy = coords[ibn].y;

        double dx = coordx - rotate_center[0];
        double dy = coordy - rotate_center[1];

        double rotCoordx = cosPsi * dx - sinPsi * dy + translation[0];
        double rotCoordy = sinPsi * dx + cosPsi * dy + translation[1];
        std::size_t matchPoint;
        double minval = std::numeric_limits<double>::max();
        for (std::size_t slave_node = slave_begn; slave_node <= slave_endn;
             ++slave_node) {
          auto slave_ibn = vertexList[slave_node].nodeIdx;

          auto slave_coordx = coords[slave_ibn].x;
          auto slave_coordy = coords[slave_ibn].y;

          double dist = 0.0;
          dist = (slave_coordx - rotCoordx) * (slave_coordx - rotCoordx) +
                 (slave_coordy - rotCoordy) * (slave_coordy - rotCoordy);
          dist = std::sqrt(dist);

          if (dist < minval) {
            minval = dist;
            matchPoint = slave_ibn;
          }
        }

        vertexList[inode].periodicPair = matchPoint;
      }
    }
    ibegn = iendn + 1;
  }

  ibegn = 0;
  for (std::size_t ib = 0; ib < numBoundSegs; ++ib) {
    std::size_t iendn = ibound[ib].bnodeIndex;
    auto name = bname[ib];
    auto type = boundaryMap.find(name)->second;
    if (type == BoundaryType::Periodic &&
        periodicMaster.find(name) != periodicMaster.end()) {
      for (std::size_t ibn = ibegn; ibn <= iendn; ++ibn) {
        std::size_t i = vertexList[ibn].nodeIdx;
        std::size_t j = vertexList[ibn].periodicPair;
        vol[i] += vol[j];
        vol[j] = vol[i];
      }
    }
    ibegn = iendn + 1;
  }
}

void Geometry::outputMeshInfo() {
  std::ofstream outputFile1("sij.txt");

  for (int i = 0; i < phyEdges; ++i) {
    outputFile1 << sij[i].x << " " << sij[i].y << " " << edge[i].nodei << " "
                << edge[i].nodej << "\n";
  }

  outputFile1.close();

  std::ofstream outputFile2("edge.txt");
  for (int i = 0; i < phyEdges; ++i) {
    outputFile2 << edge[i].nodei << " " << edge[i].nodej << "\n";
  }

  outputFile2.close();

  std::ofstream outputFile3("vol.txt");
  for (int i = 0; i < phyNodes; ++i) {
    outputFile3 << i << " " << coords[i].x << " " << coords[i].y << " "
                << vol[i] << "\n";
  }

  outputFile3.close();

  std::ofstream outputFile4("sbf.txt");
  for (int i = 0; i < numBoundFaces; ++i) {
    outputFile4 << boundaryFace[i].nodei << " " << boundaryFace[i].nodej << " "
                << sbf[i].x << " " << sbf[i].y << "\n";
  }

  outputFile4.close();

  std::ofstream outputFile5("sproj.txt");
  for (int i = 0; i < phyNodes; ++i) {
    outputFile5 << i << " " << sproj[i].x << " " << sproj[i].y << " " << vol[i]
                << "\n";
  }

  outputFile5.close();

  std::ofstream outputFile6("vertex.txt");
  for (int i = 0; i < numBoundNodes; ++i) {
    outputFile6 << i << " " << vertexList[i].nodeIdx << " "
                << vertexList[i].normal[0] << " " << vertexList[i].normal[1]
                << "\n";
  }

  outputFile6.close();

  std::ofstream outputFile7("periodic.txt");
  std::size_t ibegn = 0;
  for (std::size_t ib = 0; ib < numBoundSegs; ++ib) {
    std::size_t iendn = ibound[ib].bnodeIndex;
    auto name = bname[ib];
    auto type = boundaryMap.find(name)->second;
    if (type == BoundaryType::Periodic &&
        periodicMaster.find(name) != periodicMaster.end()) {
      for (std::size_t ibn = ibegn; ibn <= iendn; ++ibn) {
        outputFile7 << vertexList[ibn].nodeIdx << " "
                    << vertexList[ibn].periodicPair << "\n";
      }
    }
    ibegn = iendn + 1;
  }

  outputFile7.close();
}
}  // namespace preprocess
