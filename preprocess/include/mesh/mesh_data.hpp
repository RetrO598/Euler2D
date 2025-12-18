#pragma once

#include <element/element_concept.hpp>
#include <fstream>
#include <grid/point.hpp>
#include <grid/vertex.hpp>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "common/config.hpp"
#include "element/rectangle.hpp"
#include "element/triangle.hpp"
#include "element/utils.hpp"
#include "grid/edge.hpp"
#include "pre/geometry.h"
#include "pre/parameter.h"

namespace preprocess {

struct MeshDataBase {
  virtual int dim() const = 0;
  virtual ~MeshDataBase() = default;
  Index nDim;
};

template <int DIM>
struct MeshData : MeshDataBase {
  int dim() const override { return DIM; }
  Index nBound;
  std::vector<std::string> MarkerTag;
  std::vector<Point<DIM>> PointList;
  std::vector<ElementHandle> ElemList;
  std::vector<std::vector<Vertex<DIM>>> VertexList;
  std::vector<Edge<DIM>> EdgeList;
  std::vector<std::vector<ElementHandle>> BoundList;
  std::vector<std::vector<Index>> BoundElemIndex;
  std::vector<Index> nElemBound;

  void outputMeshInfo(const std::string &outFile) const;

  Geometry changeTo2Dgeometry(
      std::unordered_map<std::string, preprocess::BoundaryType> boundaryMap,
      std::vector<PeriodicInfo> periodicInfo,
      std::set<std::string> periodicMaster) const;
};

template <int DIM>
Geometry MeshData<DIM>::changeTo2Dgeometry(
    std::unordered_map<std::string, preprocess::BoundaryType> boundaryMap,
    std::vector<PeriodicInfo> periodicInfo,
    std::set<std::string> periodicMaster) const {
  Geometry geometry;
  geometry.boundaryMap = boundaryMap;
  geometry.periodicInfo = periodicInfo;
  geometry.periodicMaster = periodicMaster;
  geometry.numBoundSegs = nBound;
  geometry.bname.resize(geometry.numBoundSegs);
  geometry.ibound.resize(geometry.numBoundSegs);
  geometry.phyNodes = PointList.size();
  geometry.phyEdges = EdgeList.size();

  int numFaces = 0;
  int numNodes = 0;
  const Index nMarkers = static_cast<Index>(MarkerTag.size());
  std::vector<int> faceStart(nMarkers, 0);
  std::vector<int> nodeStart(nMarkers, 0);

  for (Index i = 0; i < nMarkers; ++i) {
    faceStart[i] = numFaces;
    nodeStart[i] = numNodes;

    numNodes += static_cast<int>(VertexList[i].size());
    numFaces += static_cast<int>(BoundList[i].size());

    geometry.ibound[i].bnodeIndex = numNodes;
    geometry.ibound[i].bfaceIndex = numFaces;
    geometry.ibound[i].bfaceIndex--;
    geometry.ibound[i].bnodeIndex--;
    geometry.bname[i] = MarkerTag[i];
  }

  geometry.numBoundFaces = numFaces;
  geometry.numBoundNodes = numNodes;

  geometry.boundaryFace.resize(numFaces);
  geometry.vertexList.resize(numNodes);

  for (Index i = 0; i < nMarkers; ++i) {
    const int fstart = faceStart[i];
    const int nstart = nodeStart[i];

    for (Index ibf = 0; ibf < BoundList[i].size(); ++ibf) {
      int idx = fstart + static_cast<int>(ibf);
      geometry.boundaryFace[idx].nodei =
          static_cast<int>(BoundList[i][ibf].getNodes(0));
      geometry.boundaryFace[idx].nodej =
          static_cast<int>(BoundList[i][ibf].getNodes(1));
    }

    for (Index ibn = 0; ibn < VertexList[i].size(); ++ibn) {
      int idx = nstart + static_cast<int>(ibn);
      geometry.vertexList[idx].nodeIdx =
          static_cast<int>(VertexList[i][ibn].getNode());
      geometry.vertexList[idx].normal[0] =
          static_cast<double>(-VertexList[i][ibn].getNormal(0));
      geometry.vertexList[idx].normal[1] =
          static_cast<double>(-VertexList[i][ibn].getNormal(1));
    }
  }

  geometry.edge.resize(EdgeList.size());
  geometry.sij.resize(EdgeList.size());
  for (Index i = 0; i < EdgeList.size(); ++i) {
    geometry.edge[i].nodei = EdgeList[i].getNode(0);
    geometry.edge[i].nodej = EdgeList[i].getNode(1);
    geometry.sij[i].x = EdgeList[i].getNormal(0);
    geometry.sij[i].y = EdgeList[i].getNormal(1);
  }

  geometry.coords.resize(PointList.size());
  geometry.vol.resize(PointList.size());
  for (Index i = 0; i < PointList.size(); ++i) {
    geometry.coords[i].x = PointList[i].getCoord(0);
    geometry.coords[i].y = PointList[i].getCoord(1);
    geometry.vol[i] = PointList[i].getVolume();
  }

  geometry.sbf.resize(geometry.numBoundFaces);
  for (Index i = 0; i < nMarkers; ++i) {
    const int fstart = faceStart[i];
    for (Index ibf = 0; ibf < BoundList[i].size(); ++ibf) {
      int idx = fstart + static_cast<int>(ibf);
      int n1 = geometry.boundaryFace[idx].nodei;
      int n2 = geometry.boundaryFace[idx].nodej;

      geometry.sbf[idx].x = geometry.coords[n2].y - geometry.coords[n1].y;
      geometry.sbf[idx].y = geometry.coords[n1].x - geometry.coords[n2].x;

      int iElem = static_cast<int>(BoundElemIndex[i][ibf]);
      if (iElem < 0 || static_cast<size_t>(iElem) >= ElemList.size()) {
        continue;
      }

      double cx = ElemList[iElem].getCoordCenter(0);
      double cy = ElemList[iElem].getCoordCenter(1);
      double xm = 0.5 * (geometry.coords[n1].x + geometry.coords[n2].x);
      double ym = 0.5 * (geometry.coords[n1].y + geometry.coords[n2].y);
      double vprod =
          geometry.sbf[idx].x * (xm - cx) + geometry.sbf[idx].y * (ym - cy);

      if (vprod < 0.0) {
        geometry.sbf[idx].x = -geometry.sbf[idx].x;
        geometry.sbf[idx].y = -geometry.sbf[idx].y;
      }
    }
  }

  geometry.elements.reserve(ElemList.size());
  geometry.numElems = ElemList.size();
  for (Index i = 0; i < ElemList.size(); ++i) {
    if (ElemList[i].getVTKtype() == TRIANGLE) {
      geometry.elements.push_back(Triangle<2>(ElemList[i].getNodes(0),
                                              ElemList[i].getNodes(1),
                                              ElemList[i].getNodes(2)));
    } else if (ElemList[i].getVTKtype() == RECTANGLE) {
      geometry.elements.push_back(
          Rectangle<2>(ElemList[i].getNodes(0), ElemList[i].getNodes(1),
                       ElemList[i].getNodes(2), ElemList[i].getNodes(3)));
    }
  }

  geometry.pointList = PointList;

  geometry.sproj.resize(geometry.phyNodes);
  geometry.MatchPeriodic();
  geometry.CheckMetrics();
  geometry.FaceVectorsSymm();
  geometry.volumeProjections();

  return geometry;
}

template <int DIM>
void MeshData<DIM>::outputMeshInfo(const std::string &outFile) const {
  std::ofstream out(outFile);
  if (!out.good()) {
    return;
  }
  out << "================== POINTS INFORMATION ==================\n";
  out << "# nPoint nPointDomain nDim\n";

  out << PointList.size() << " " << PointList.size() << " " << nDim << "\n";
  out << "# PointIndex Domain Coordinates... NeighborElems NeighborPoints "
         "EdgeIndices\n";
  for (Index i = 0; i < PointList.size(); ++i) {
    out << i << " " << "1" << " ";
    for (int d = 0; d < nDim; ++d) {
      out << PointList[i].getCoord(d) << " ";
    }
    for (Index j = 0; j < PointList[i].getnElem(); ++j) {
      out << PointList[i].getElem(j) << " ";
    }
    for (Index k = 0; k < PointList[i].getnPoint(); ++k) {
      out << PointList[i].getPoint(k) << " ";
    }
    for (Index k = 0; k < PointList[i].getnPoint(); ++k) {
      auto iedge = PointList[i].getEdge(k);
      if (iedge == INVALID_INDEX) {
        out << "-1 ";
      } else {
        out << iedge << " ";
      }
    }
    out << "\n";
  }

  out << "\n================== ELEMENTS INFORMATION ==================\n";
  out << "# nElem nDim\n";
  out << ElemList.size() << " " << nDim << "\n";
  out << "# ElementIndex VTKType NodeIndices...\n";
  for (Index i = 0; i < ElemList.size(); ++i) {
    out << i << " ";
    auto vtk = ElemList[i].getVTKtype();
    out << vtk << " ";
    Index nNodes = ElemList[i].getnNodes();
    for (Index j = 0; j < nNodes; ++j) {
      out << ElemList[i].getNodes(j) << " ";
    }
    out << "\n";
  }

  out << "\n================== CONTROL VOLUMES INFORMATION "
         "==================\n";
  out << "# PointIndex nSurroundingPoints SurroundingPoints... Volume\n";
  for (Index i = 0; i < PointList.size(); ++i) {
    out << i << " ";
    Index nSurroundingPoints = PointList[i].getnPoint();
    out << nSurroundingPoints << " ";
    for (Index j = 0; j < nSurroundingPoints; ++j) {
      out << PointList[i].getPoint(j) << " ";
    }
    out << PointList[i].getVolume() << "\n";
  }

  out << "\n================== EDGE INFORMATION "
         "==================\n";
  out << "# nEdge nDim\n";
  out << EdgeList.size() << " " << nDim << "\n";
  out << "# EdgeIndex Node0 Node1 EdgeCenterCoords... Normal...\n";
  for (Index i = 0; i < EdgeList.size(); ++i) {
    out << i << " ";
    out << EdgeList[i].getNode(0) << " " << EdgeList[i].getNode(1) << " ";
    for (Index d = 0; d < nDim; ++d) {
      out << EdgeList[i].getCoord(d) << " ";
    }
    auto normal = EdgeList[i].getNormal();
    for (Index d = 0; d < nDim; ++d) {
      out << normal[d] << " ";
    }
    out << "\n";
  }

  out << "\n================== BOUNDARY INFORMATION "
         "==================\n";
  for (Index i = 0; i < nBound; ++i) {
    out << "# MarkerIndex MarkerTag nElemBound\n";
    out << i << " " << MarkerTag[i] << " " << nElemBound[i] << "\n";
    out << "# BoundElementIndex VTKType NodeIndices...\n";
    for (Index j = 0; j < nElemBound[i]; ++j) {
      out << j << " ";
      auto vtk = BoundList[i][j].getVTKtype();
      out << vtk << " ";
      Index nNodes = BoundList[i][j].getnNodes();
      for (Index k = 0; k < nNodes; ++k) {
        out << BoundList[i][j].getNodes(k) << " ";
      }
      out << "\n";
    }
  }

  out << "\n================== VERTEX INFORMATION "
         "==================\n";
  if (!VertexList.empty()) {
    for (Index i = 0; i < nBound; ++i) {
      auto nVertex = static_cast<Index>(VertexList[i].size());
      out << "# MarkerTag nVertex\n";
      out << MarkerTag[i] << " " << nVertex << "\n";
      out << "# VertexIndex PointIndex VertexCoord... Normal... "
             "NormalNeighbor\n";
      for (Index j = 0; j < nVertex; ++j) {
        out << j << " " << VertexList[i][j].getNode() << " ";
        for (Index d = 0; d < nDim; ++d) {
          out << VertexList[i][j].getCoord(d) << " ";
        }
        auto normal = VertexList[i][j].getNormal();
        for (Index d = 0; d < nDim; ++d) {
          out << normal[d] << " ";
        }
        if (VertexList[i][j].getNormalNeighbor() != INVALID_INDEX) {
          out << VertexList[i][j].getNormalNeighbor() << "\n";
        } else {
          out << "0" << "\n";
        }
      }
    }
  }
  out.close();
};
}  // namespace preprocess