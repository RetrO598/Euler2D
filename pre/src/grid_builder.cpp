#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>

#include "common/config.hpp"
#include "element/utils.hpp"
#include "grid/grid_builder.hpp"
#include "mesh/mesh_data.hpp"
namespace preprocess {
namespace GridBuilder {

template <int DIM>
void FVMBuilder::buildImpl(MeshData<DIM> &mesh) {
  if constexpr (DIM == 2) {
    std::cout << "Building 2D grid from " << mesh.PointList.size() << " points."
              << std::endl;
    SetPointConnectivity(mesh);

    // for (Index iPoint = 0; iPoint < mesh.PointList.size(); ++iPoint) {
    //   mesh.PointList[iPoint].resetElem();
    //   mesh.PointList[iPoint].resetPoint();
    //   mesh.PointList[iPoint].resetBoundary();
    // }
    // auto &boundList = mesh.BoundList;
    // for (Index iMarker = 0; iMarker < mesh.nBound; ++iMarker) {
    //   for (Index iElem = 0; iElem < mesh.nElemBound[iMarker]; ++iElem) {
    //     std::string markerTag = mesh.MarkerTag[iMarker];
    //     for (Index iNode = 0; iNode < boundList[iMarker][iElem].getnNodes();
    //          ++iNode) {
    //       Index iPoint = boundList[iMarker][iElem].getNodes(iNode);
    //       mesh.PointList[iPoint].setBoundary(mesh.nBound);
    //       mesh.PointList[iPoint].setPhysicalBoundary(true);
    //       // if (markerTag == "EULER_WALL" || markerTag == "NOSLIP_WALL") {
    //       //   pointList[invResult[iPoint]].setSolidBoundary(true);
    //       // }
    //     }
    //   }
    // }
    SetRCMOrdering(mesh);
    SetPointConnectivity(mesh);
    SetElementConnectivity(mesh);
    SetBoundVolume(mesh);
    CheckOrientation(mesh);
    SetEdges(mesh);
    SetVertex(mesh);
    SetCG(mesh);
    SetControlVolume(mesh);
    SetBoundControlVolume(mesh);
    FindNormalNeighbor(mesh);
  } else if constexpr (DIM == 3) {
    std::cout << "Building 3D grid from " << mesh.PointList.size() << " points."
              << std::endl;
    SetPointConnectivity(mesh);
    SetRCMOrdering(mesh);
    SetPointConnectivity(mesh);
    SetElementConnectivity(mesh);
    SetBoundVolume(mesh);
    CheckOrientation(mesh);
    SetEdges(mesh);
    SetVertex(mesh);
    SetCG(mesh);
    SetControlVolume(mesh);
    SetBoundControlVolume(mesh);
    FindNormalNeighbor(mesh);
  }
}

void FVMBuilder::build(MeshDataBase &mesh) {
  if (mesh.dim() == 2) {
    auto &mesh2d = dynamic_cast<MeshData<2> &>(mesh);
    buildImpl(mesh2d);
  } else if (mesh.dim() == 3) {
    auto &mesh3d = dynamic_cast<MeshData<3> &>(mesh);
    buildImpl(mesh3d);
  } else {
    throw std::runtime_error("Unsupported mesh dimension in FVMBuilder.");
  }
}

template <int DIM>
void FVMBuilder::SetPointConnectivity(MeshData<DIM> &mesh) {
  auto &pointList = mesh.PointList;
  auto &elemList = mesh.ElemList;
  Index Point_Neighbor;

  for (Index iElem = 0; iElem < elemList.size(); ++iElem) {
    for (Index iNode = 0; iNode < elemList[iElem].getnNodes(); ++iNode) {
      auto iPoint = elemList[iElem].getNodes(iNode);
      pointList[iPoint].setElem(iElem);
    }
  }

  for (Index iPoint = 0; iPoint < pointList.size(); ++iPoint) {
    for (Index iElem = 0; iElem < pointList[iPoint].getnElem(); ++iElem) {
      auto jElem = pointList[iPoint].getElem(iElem);

      for (Index iNode = 0; iNode < elemList[jElem].getnNodes(); ++iNode) {
        if (elemList[jElem].getNodes(iNode) == iPoint) {
          // Use getnNeighborNodes(iNode) to get the count of neighboring nodes
          // for this specific node
          for (Index iNeighbor = 0;
               iNeighbor < elemList[jElem].getnNeighborNodes(iNode);
               ++iNeighbor) {
            auto neighborNodeIdx =
                elemList[jElem].getNodesNeighbor(iNode, iNeighbor);
            Point_Neighbor = elemList[jElem].getNodes(neighborNodeIdx);
            pointList[iPoint].setPoint(Point_Neighbor);
          }
        }
      }
    }
  }
}

template <int DIM>
void FVMBuilder::SetRCMOrdering(MeshData<DIM> &mesh) {
  Index nPoints = mesh.PointList.size();
  std::vector<Index> result;
  std::vector<Index> auxQueue;
  std::vector<Index> Queue;
  auto &pointList = mesh.PointList;
  auto &elemList = mesh.ElemList;
  std::unique_ptr<bool[]> inQueue = std::make_unique<bool[]>(nPoints);

  for (Index iPoint = 0; iPoint < nPoints; ++iPoint) {
    inQueue[iPoint] = false;
  }

  Index degree = 0;
  Index minDegree = pointList[0].getnPoint();
  Index addPoint = 0;
  for (Index iPoint = 1; iPoint < nPoints; ++iPoint) {
    degree = pointList[iPoint].getnPoint();
    if (degree < minDegree) {
      minDegree = degree;
      addPoint = iPoint;
    }
  }

  result.push_back(addPoint);
  inQueue[addPoint] = true;

  do {
    auxQueue.clear();
    for (Index iNode = 0; iNode < pointList[addPoint].getnPoint(); ++iNode) {
      Index adjPoint = pointList[addPoint].getPoint(iNode);
      if ((!inQueue[adjPoint]) && (adjPoint < nPoints)) {
        auxQueue.push_back(adjPoint);
      }
    }

    if (auxQueue.size() != 0) {
      for (Index iNode = 0; iNode < auxQueue.size(); ++iNode) {
        for (Index jNode = 0; jNode < auxQueue.size() - 1 - iNode; ++jNode) {
          if (pointList[auxQueue[jNode]].getnPoint() >
              pointList[auxQueue[jNode + 1]].getnPoint()) {
            Index temp = auxQueue[jNode];
            auxQueue[jNode] = auxQueue[jNode + 1];
            auxQueue[jNode + 1] = temp;
          }
        }
      }

      Queue.insert(Queue.end(), auxQueue.begin(), auxQueue.end());
      for (Index iNode = 0; iNode < auxQueue.size(); ++iNode) {
        inQueue[auxQueue[iNode]] = true;
      }
    }

    if (Queue.size() != 0) {
      addPoint = Queue[0];
      result.push_back(Queue[0]);
      Queue.erase(Queue.begin(), Queue.begin() + 1);
    }

  } while (Queue.size() != 0);

  for (Index iPoint = 0; iPoint < nPoints; ++iPoint) {
    if (inQueue[iPoint] == false) {
      result.push_back(iPoint);
    }
  }

  std::reverse(result.begin(), result.end());

  for (Index iPoint = 0; iPoint < nPoints; ++iPoint) {
    pointList[iPoint].resetElem();
    pointList[iPoint].resetPoint();
    pointList[iPoint].resetBoundary();
    pointList[iPoint].setPhysicalBoundary(false);
    pointList[iPoint].setSolidBoundary(false);
    pointList[iPoint].setDomain(true);
  }

  std::vector<Index> auxPointIndex(nPoints);
  std::vector<std::array<real, DIM>> auxCoord(nPoints);
  for (Index iPoint = 0; iPoint < nPoints; ++iPoint) {
    auxPointIndex[iPoint] = pointList[iPoint].getIndex();
    auxCoord[iPoint] = pointList[iPoint].getCoord();
  }

  for (Index iPoint = 0; iPoint < nPoints; ++iPoint) {
    pointList[iPoint].setIndex(auxPointIndex[result[iPoint]]);
    pointList[iPoint].setCoord(auxCoord[result[iPoint]]);
  }

  std::vector<Index> invResult(nPoints);
  for (Index iPoint = 0; iPoint < nPoints; ++iPoint) {
    invResult[result[iPoint]] = iPoint;
  }

  for (Index iElem = 0; iElem < elemList.size(); ++iElem) {
    for (Index iNode = 0; iNode < elemList[iElem].getnNodes(); ++iNode) {
      Index iPoint = elemList[iElem].getNodes(iNode);
      elemList[iElem].setNodes(iNode, invResult[iPoint]);
    }
  }

  auto &boundList = mesh.BoundList;
  for (Index iMarker = 0; iMarker < mesh.nBound; ++iMarker) {
    for (Index iElem = 0; iElem < mesh.nElemBound[iMarker]; ++iElem) {
      std::string markerTag = mesh.MarkerTag[iMarker];
      for (Index iNode = 0; iNode < boundList[iMarker][iElem].getnNodes();
           ++iNode) {
        Index iPoint = boundList[iMarker][iElem].getNodes(iNode);
        boundList[iMarker][iElem].setNodes(iNode, invResult[iPoint]);
        pointList[invResult[iPoint]].setBoundary(mesh.nBound);
        pointList[invResult[iPoint]].setPhysicalBoundary(true);
        if (markerTag == "EULER_WALL" || markerTag == "NOSLIP_WALL") {
          pointList[invResult[iPoint]].setSolidBoundary(true);
        }
      }
    }
  }
}

template <int DIM>
void FVMBuilder::SetElementConnectivity(MeshData<DIM> &mesh) {
  auto &elemList = mesh.ElemList;
  auto &pointList = mesh.PointList;
  Index nElem = elemList.size();
  Index firstElemFace, secondElemFace;

  for (Index iElem = 0; iElem < nElem; ++iElem) {
    for (Index iFace = 0; iFace < elemList[iElem].getnFaces(); ++iFace) {
      for (Index iNode = 0; iNode < elemList[iElem].getnNodesFace(iFace);
           ++iNode) {
        Index facePoint =
            elemList[iElem].getNodes(elemList[iElem].getFaces(iFace, iNode));

        for (Index jElem = 0; jElem < pointList[facePoint].getnElem();
             ++jElem) {
          Index testElem = pointList[facePoint].getElem(jElem);
          if ((elemList[iElem].getNeighborElems(iFace) == INVALID_INDEX) &&
              (iElem < testElem) &&
              (findFace(mesh, iElem, testElem, firstElemFace,
                        secondElemFace))) {
            elemList[iElem].setNeighborElems(firstElemFace, testElem);
            elemList[testElem].setNeighborElems(secondElemFace, iElem);
          }
        }
      }
    }
  }
}

template <int DIM>
bool FVMBuilder::findFace(MeshData<DIM> &mesh, Index elem1, Index elem2,
                          Index &elem1Face, Index &elem2Face) {
  if (elem1 == elem2) {
    return false;
  }

  auto &elemList = mesh.ElemList;
  std::vector<Index> commonPoints;
  std::vector<Index> pointFaceFirst;
  std::vector<Index> pointFaceSecond;
  Index iPoint{}, jPoint{};
  bool faceFirstFound = false;
  bool faceSecondFound = false;

  Index kNode = 0;
  for (Index iNode = 0; iNode < elemList[elem1].getnNodes(); ++iNode) {
    iPoint = elemList[elem1].getNodes(iNode);
    for (Index jNode = 0; jNode < elemList[elem2].getnNodes(); ++jNode) {
      jPoint = elemList[elem2].getNodes(jNode);
      if (iPoint == jPoint) {
        commonPoints.push_back(iPoint);
        break;
      }
    }
  }

  std::sort(commonPoints.begin(), commonPoints.end());
  auto iterPoint = std::unique(commonPoints.begin(), commonPoints.end());
  commonPoints.resize(std::distance(commonPoints.begin(), iterPoint));

  for (Index iFace = 0; iFace < elemList[elem1].getnFaces(); ++iFace) {
    Index nNodesFace = elemList[elem1].getnNodesFace(iFace);
    for (Index iNode = 0; iNode < nNodesFace; ++iNode) {
      Index faceNode = elemList[elem1].getFaces(iFace, iNode);
      pointFaceFirst.push_back(elemList[elem1].getNodes(faceNode));
    }

    std::sort(pointFaceFirst.begin(), pointFaceFirst.end());

    auto mypair = std::mismatch(pointFaceFirst.begin(), pointFaceFirst.end(),
                                commonPoints.begin());
    if (mypair.first == pointFaceFirst.end()) {
      elem1Face = iFace;
      faceFirstFound = true;
      break;
    }

    pointFaceFirst.erase(pointFaceFirst.begin(), pointFaceFirst.end());
  }

  for (Index iFace = 0; iFace < elemList[elem2].getnFaces(); ++iFace) {
    Index nNodesFace = elemList[elem2].getnNodesFace(iFace);
    for (Index iNode = 0; iNode < nNodesFace; ++iNode) {
      Index faceNode = elemList[elem2].getFaces(iFace, iNode);
      pointFaceSecond.push_back(elemList[elem2].getNodes(faceNode));
    }

    std::sort(pointFaceSecond.begin(), pointFaceSecond.end());

    auto mypair = std::mismatch(pointFaceSecond.begin(), pointFaceSecond.end(),
                                commonPoints.begin());
    if (mypair.first == pointFaceSecond.end()) {
      elem2Face = iFace;
      faceSecondFound = true;
      break;
    }

    pointFaceSecond.erase(pointFaceSecond.begin(), pointFaceSecond.end());
  }

  if (faceFirstFound && faceSecondFound) {
    return true;
  } else {
    return false;
  }
}

template <int DIM>
void FVMBuilder::SetBoundVolume(MeshData<DIM> &mesh) {
  auto &pointList = mesh.PointList;
  auto &boundList = mesh.BoundList;
  auto &elemList = mesh.ElemList;
  bool checkVol;

  for (Index iMarker = 0; iMarker < mesh.nBound; ++iMarker) {
    for (Index iElemSurface = 0; iElemSurface < mesh.nElemBound[iMarker];
         ++iElemSurface) {
      Index point = boundList[iMarker][iElemSurface].getNodes(0);
      checkVol = false;

      for (Index iElem = 0; iElem < pointList[point].getnElem(); ++iElem) {
        Index count = 0;
        Index iElemDomain = pointList[point].getElem(iElem);
        for (Index iNodeDomain = 0;
             iNodeDomain < elemList[iElemDomain].getnNodes(); ++iNodeDomain) {
          Index pointDomain = elemList[iElemDomain].getNodes(iNodeDomain);
          for (Index iNodeSurface = 0;
               iNodeSurface < boundList[iMarker][iElemSurface].getnNodes();
               ++iNodeSurface) {
            Index pointSurface =
                boundList[iMarker][iElemSurface].getNodes(iNodeSurface);
            if (pointDomain == pointSurface) {
              count++;
            }
            if (count == boundList[iMarker][iElemSurface].getnNodes()) {
              break;
            }
          }
          if (count == boundList[iMarker][iElemSurface].getnNodes()) {
            break;
          }
        }
        if (count == boundList[iMarker][iElemSurface].getnNodes()) {
          mesh.BoundElemIndex[iMarker][iElemSurface] = iElemDomain;
          checkVol = true;
          break;
        }
      }
      if (!checkVol) {
        std::cout << "The surface element (" << iMarker << ", " << iElemSurface
                  << ") doesn't have an associated volume element."
                  << std::endl;
        std::exit(1);
      }
    }
  }
}

template <int DIM>
void FVMBuilder::CheckOrientation(MeshData<DIM> &mesh) {
  Index nElem = mesh.ElemList.size();
  auto &pointList = mesh.PointList;
  auto &elemList = mesh.ElemList;
  Index point1, point2, point3, point4, point5, point6;
  Index point1Surface, point2Surface, point3Surface, point4Surface;
  std::array<real, DIM> coord1, coord2, coord3, coord4, coord5, coord6, a, b, c,
      n;
  real test1, test2, test3, test4;

  for (Index iElem = 0; iElem < nElem; ++iElem) {
    if (elemList[iElem].getVTKtype() == TRIANGLE) {
      point1 = elemList[iElem].getNodes(0);
      coord1 = mesh.PointList[point1].getCoord();
      point2 = elemList[iElem].getNodes(1);
      coord2 = mesh.PointList[point2].getCoord();
      point3 = elemList[iElem].getNodes(2);
      coord3 = mesh.PointList[point3].getCoord();
      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord2[i] - coord1[i]);
        b[i] = 0.5 * (coord3[i] - coord1[i]);
      }
      real test = a[0] * b[1] - b[0] * a[1];
      if (test < 0.0) {
        elemList[iElem].changeOrientation();
      }
    }

    if (elemList[iElem].getVTKtype() == RECTANGLE) {
      point1 = elemList[iElem].getNodes(0);
      coord1 = mesh.PointList[point1].getCoord();
      point2 = elemList[iElem].getNodes(1);
      coord2 = mesh.PointList[point2].getCoord();
      point3 = elemList[iElem].getNodes(2);
      coord3 = mesh.PointList[point3].getCoord();
      point4 = elemList[iElem].getNodes(3);
      coord4 = mesh.PointList[point4].getCoord();

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord2[i] - coord1[i]);
        b[i] = 0.5 * (coord3[i] - coord1[i]);
      }
      test1 = a[0] * b[1] - b[0] * a[1];

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord3[i] - coord2[i]);
        b[i] = 0.5 * (coord4[i] - coord2[i]);
      }
      test2 = a[0] * b[1] - b[0] * a[1];

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord4[i] - coord3[i]);
        b[i] = 0.5 * (coord1[i] - coord3[i]);
      }
      test3 = a[0] * b[1] - b[0] * a[1];

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord1[i] - coord4[i]);
        b[i] = 0.5 * (coord3[i] - coord4[i]);
      }
      test4 = a[0] * b[1] - b[0] * a[1];

      if ((test1 < 0.0) && (test2 < 0.0) && (test3 < 0.0) && (test4 < 0.0))
        elemList[iElem].changeOrientation();
    }

    if (elemList[iElem].getVTKtype() == TETRAHEDRON) {
      point1 = elemList[iElem].getNodes(0);
      coord1 = mesh.PointList[point1].getCoord();
      point2 = elemList[iElem].getNodes(1);
      coord2 = mesh.PointList[point2].getCoord();
      point3 = elemList[iElem].getNodes(2);
      coord3 = mesh.PointList[point3].getCoord();
      point4 = elemList[iElem].getNodes(3);
      coord4 = mesh.PointList[point4].getCoord();

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord2[i] - coord1[i]);
        b[i] = 0.5 * (coord3[i] - coord1[i]);
        c[i] = coord4[i] - coord1[i];
      }
      n[0] = a[1] * b[2] - b[1] * a[2];
      n[1] = -(a[0] * b[2] - b[0] * a[2]);
      n[2] = a[0] * b[1] - b[0] * a[1];

      real test = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];
      if (test < 0.0) elemList[iElem].changeOrientation();
    }

    if (elemList[iElem].getVTKtype() == WEDGE) {
      point1 = elemList[iElem].getNodes(0);
      coord1 = mesh.PointList[point1].getCoord();
      point2 = elemList[iElem].getNodes(1);
      coord2 = mesh.PointList[point2].getCoord();
      point3 = elemList[iElem].getNodes(2);
      coord3 = mesh.PointList[point3].getCoord();
      point4 = elemList[iElem].getNodes(3);
      coord4 = mesh.PointList[point4].getCoord();
      point5 = elemList[iElem].getNodes(4);
      coord5 = mesh.PointList[point5].getCoord();
      point6 = elemList[iElem].getNodes(5);
      coord6 = mesh.PointList[point6].getCoord();

      for (Index iDim = 0; iDim < DIM; iDim++) {
        a[iDim] = 0.5 * (coord3[iDim] - coord1[iDim]);
        b[iDim] = 0.5 * (coord2[iDim] - coord1[iDim]);
        c[iDim] = (coord4[iDim] - coord1[iDim]) +
                  (coord5[iDim] - coord2[iDim]) + (coord6[iDim] - coord3[iDim]);
      }

      /*--- The normal vector should point to the interior of the element ---*/
      n[0] = a[1] * b[2] - b[1] * a[2];
      n[1] = -(a[0] * b[2] - b[0] * a[2]);
      n[2] = a[0] * b[1] - b[0] * a[1];

      test1 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

      for (Index iDim = 0; iDim < DIM; iDim++) {
        a[iDim] = 0.5 * (coord5[iDim] - coord4[iDim]);
        b[iDim] = 0.5 * (coord6[iDim] - coord4[iDim]);
        c[iDim] = (coord1[iDim] - coord4[iDim]) +
                  (coord2[iDim] - coord5[iDim]) + (coord3[iDim] - coord6[iDim]);
      }

      /*--- The normal vector should point to the interior of the element ---*/
      n[0] = a[1] * b[2] - b[1] * a[2];
      n[1] = -(a[0] * b[2] - b[0] * a[2]);
      n[2] = a[0] * b[1] - b[0] * a[1];

      test2 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

      if ((test1 < 0.0) || (test2 < 0.0)) elemList[iElem].changeOrientation();
    }

    if (elemList[iElem].getVTKtype() == HEXAHEDRON) {
      point1 = elemList[iElem].getNodes(0);
      coord1 = pointList[point1].getCoord();
      point2 = elemList[iElem].getNodes(1);
      coord2 = pointList[point2].getCoord();
      point3 = elemList[iElem].getNodes(2);
      coord3 = pointList[point3].getCoord();
      point4 = elemList[iElem].getNodes(5);
      coord4 = pointList[point4].getCoord();

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord2[i] - coord1[i]);
        b[i] = 0.5 * (coord3[i] - coord1[i]);
        c[i] = coord4[i] - coord1[i];
      }
      n[0] = a[1] * b[2] - b[1] * a[2];
      n[1] = -(a[0] * b[2] - b[0] * a[2]);
      n[2] = a[0] * b[1] - b[0] * a[1];

      test1 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

      point1 = elemList[iElem].getNodes(2);
      coord1 = pointList[point1].getCoord();
      point2 = elemList[iElem].getNodes(3);
      coord2 = pointList[point2].getCoord();
      point3 = elemList[iElem].getNodes(0);
      coord3 = pointList[point3].getCoord();
      point4 = elemList[iElem].getNodes(7);
      coord4 = pointList[point4].getCoord();

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord2[i] - coord1[i]);
        b[i] = 0.5 * (coord3[i] - coord1[i]);
        c[i] = coord4[i] - coord1[i];
      }
      n[0] = a[1] * b[2] - b[1] * a[2];
      n[1] = -(a[0] * b[2] - b[0] * a[2]);
      n[2] = a[0] * b[1] - b[0] * a[1];

      test2 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

      point1 = elemList[iElem].getNodes(1);
      coord1 = pointList[point1].getCoord();
      point2 = elemList[iElem].getNodes(2);
      coord2 = pointList[point2].getCoord();
      point3 = elemList[iElem].getNodes(3);
      coord3 = pointList[point3].getCoord();
      point4 = elemList[iElem].getNodes(6);
      coord4 = pointList[point4].getCoord();

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord2[i] - coord1[i]);
        b[i] = 0.5 * (coord3[i] - coord1[i]);
        c[i] = coord4[i] - coord1[i];
      }
      n[0] = a[1] * b[2] - b[1] * a[2];
      n[1] = -(a[0] * b[2] - b[0] * a[2]);
      n[2] = a[0] * b[1] - b[0] * a[1];

      test3 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

      point1 = elemList[iElem].getNodes(3);
      coord1 = pointList[point1].getCoord();
      point2 = elemList[iElem].getNodes(0);
      coord2 = pointList[point2].getCoord();
      point3 = elemList[iElem].getNodes(1);
      coord3 = pointList[point3].getCoord();
      point4 = elemList[iElem].getNodes(4);
      coord4 = pointList[point4].getCoord();

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord2[i] - coord1[i]);
        b[i] = 0.5 * (coord3[i] - coord1[i]);
        c[i] = coord4[i] - coord1[i];
      }
      n[0] = a[1] * b[2] - b[1] * a[2];
      n[1] = -(a[0] * b[2] - b[0] * a[2]);
      n[2] = a[0] * b[1] - b[0] * a[1];

      test4 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

      if ((test1 < 0.0) || (test2 < 0.0) || (test3 < 0.0) || (test4 < 0.0))
        elemList[iElem].changeOrientation();
    }

    if (elemList[iElem].getVTKtype() == PYRAMID) {
      point1 = elemList[iElem].getNodes(0);
      coord1 = pointList[point1].getCoord();
      point2 = elemList[iElem].getNodes(1);
      coord2 = pointList[point2].getCoord();
      point3 = elemList[iElem].getNodes(2);
      coord3 = pointList[point3].getCoord();
      point4 = elemList[iElem].getNodes(4);
      coord4 = pointList[point4].getCoord();

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord2[i] - coord1[i]);
        b[i] = 0.5 * (coord3[i] - coord1[i]);
        c[i] = coord4[i] - coord1[i];
      }
      n[0] = a[1] * b[2] - b[1] * a[2];
      n[1] = -(a[0] * b[2] - b[0] * a[2]);
      n[2] = a[0] * b[1] - b[0] * a[1];

      test1 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

      point1 = elemList[iElem].getNodes(2);
      coord1 = pointList[point1].getCoord();
      point2 = elemList[iElem].getNodes(3);
      coord2 = pointList[point2].getCoord();
      point3 = elemList[iElem].getNodes(0);
      coord3 = pointList[point3].getCoord();
      point4 = elemList[iElem].getNodes(4);
      coord4 = pointList[point4].getCoord();

      for (Index i = 0; i < DIM; ++i) {
        a[i] = 0.5 * (coord2[i] - coord1[i]);
        b[i] = 0.5 * (coord3[i] - coord1[i]);
        c[i] = coord4[i] - coord1[i];
      }
      n[0] = a[1] * b[2] - b[1] * a[2];
      n[1] = -(a[0] * b[2] - b[0] * a[2]);
      n[2] = a[0] * b[1] - b[0] * a[1];

      test2 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

      if ((test1 < 0.0) || (test2 < 0.0)) elemList[iElem].changeOrientation();
    }
  }

  bool find;
  Index iElemDomain, pointDomain, pointSurface;
  for (Index iMarker = 0; iMarker < mesh.nBound; ++iMarker) {
    for (Index iElemSurface = 0; iElemSurface < mesh.nElemBound[iMarker];
         ++iElemSurface) {
      iElemDomain = mesh.BoundElemIndex[iMarker][iElemSurface];
      for (Index iNodeDomain = 0;
           iNodeDomain < elemList[iElemDomain].getnNodes(); ++iNodeDomain) {
        pointDomain = elemList[iElemDomain].getNodes(iNodeDomain);
        find = false;
        for (Index iNodeSurface = 0;
             iNodeSurface < mesh.BoundList[iMarker][iElemSurface].getnNodes();
             ++iNodeSurface) {
          pointSurface =
              mesh.BoundList[iMarker][iElemSurface].getNodes(iNodeSurface);
          if (pointSurface == pointDomain) {
            find = true;
            break;
          }
        }
        if (!find) {
          break;
        }
      }

      if (mesh.BoundList[iMarker][iElemSurface].getVTKtype() == LINE) {
        point1Surface = mesh.BoundList[iMarker][iElemSurface].getNodes(0);
        coord1 = pointList[point1Surface].getCoord();
        point2Surface = mesh.BoundList[iMarker][iElemSurface].getNodes(1);
        coord2 = pointList[point2Surface].getCoord();
        coord3 = pointList[pointDomain].getCoord();

        for (Index i = 0; i < DIM; i++) {
          a[i] = 0.5 * (coord2[i] - coord1[i]);
          b[i] = 0.5 * (coord3[i] - coord1[i]);
        }
        real test = a[0] * b[1] - b[0] * a[1];

        if (test < 0.0)
          mesh.BoundList[iMarker][iElemSurface].changeOrientation();
      }

      if (mesh.BoundList[iMarker][iElemSurface].getVTKtype() == TRIANGLE) {
        point1Surface = mesh.BoundList[iMarker][iElemSurface].getNodes(0);
        coord1 = pointList[point1Surface].getCoord();
        point2Surface = mesh.BoundList[iMarker][iElemSurface].getNodes(1);
        coord2 = pointList[point2Surface].getCoord();
        point3Surface = mesh.BoundList[iMarker][iElemSurface].getNodes(2);
        coord3 = pointList[point3Surface].getCoord();
        coord4 = pointList[pointDomain].getCoord();

        for (Index iDim = 0; iDim < DIM; iDim++) {
          a[iDim] = 0.5 * (coord2[iDim] - coord1[iDim]);
          b[iDim] = 0.5 * (coord3[iDim] - coord1[iDim]);
          c[iDim] = coord4[iDim] - coord1[iDim];
        }
        n[0] = a[1] * b[2] - b[1] * a[2];
        n[1] = -(a[0] * b[2] - b[0] * a[2]);
        n[2] = a[0] * b[1] - b[0] * a[1];

        real test = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];
        if (test < 0.0)
          mesh.BoundList[iMarker][iElemSurface].changeOrientation();
      }

      if (mesh.BoundList[iMarker][iElemSurface].getVTKtype() == RECTANGLE) {
        point1Surface = mesh.BoundList[iMarker][iElemSurface].getNodes(0);
        coord1 = pointList[point1Surface].getCoord();
        point2Surface = mesh.BoundList[iMarker][iElemSurface].getNodes(1);
        coord2 = pointList[point2Surface].getCoord();
        point3Surface = mesh.BoundList[iMarker][iElemSurface].getNodes(2);
        coord3 = pointList[point3Surface].getCoord();
        point4Surface = mesh.BoundList[iMarker][iElemSurface].getNodes(3);
        coord4 = pointList[point4Surface].getCoord();
        coord5 = pointList[pointDomain].getCoord();

        for (Index iDim = 0; iDim < DIM; iDim++) {
          a[iDim] = 0.5 * (coord2[iDim] - coord1[iDim]);
          b[iDim] = 0.5 * (coord3[iDim] - coord1[iDim]);
          c[iDim] = coord5[iDim] - coord1[iDim];
        }
        n[0] = a[1] * b[2] - b[1] * a[2];
        n[1] = -(a[0] * b[2] - b[0] * a[2]);
        n[2] = a[0] * b[1] - b[0] * a[1];
        test1 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

        for (Index iDim = 0; iDim < DIM; iDim++) {
          a[iDim] = 0.5 * (coord3[iDim] - coord2[iDim]);
          b[iDim] = 0.5 * (coord4[iDim] - coord2[iDim]);
          c[iDim] = coord5[iDim] - coord2[iDim];
        }
        n[0] = a[1] * b[2] - b[1] * a[2];
        n[1] = -(a[0] * b[2] - b[0] * a[2]);
        n[2] = a[0] * b[1] - b[0] * a[1];
        test2 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

        for (Index iDim = 0; iDim < DIM; iDim++) {
          a[iDim] = 0.5 * (coord4[iDim] - coord3[iDim]);
          b[iDim] = 0.5 * (coord1[iDim] - coord3[iDim]);
          c[iDim] = coord5[iDim] - coord3[iDim];
        }
        n[0] = a[1] * b[2] - b[1] * a[2];
        n[1] = -(a[0] * b[2] - b[0] * a[2]);
        n[2] = a[0] * b[1] - b[0] * a[1];
        test3 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

        for (Index iDim = 0; iDim < DIM; iDim++) {
          a[iDim] = 0.5 * (coord1[iDim] - coord4[iDim]);
          b[iDim] = 0.5 * (coord3[iDim] - coord4[iDim]);
          c[iDim] = coord5[iDim] - coord4[iDim];
        }
        n[0] = a[1] * b[2] - b[1] * a[2];
        n[1] = -(a[0] * b[2] - b[0] * a[2]);
        n[2] = a[0] * b[1] - b[0] * a[1];
        test4 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

        if ((test1 < 0.0) && (test2 < 0.0) && (test3 < 0.0) && (test4 < 0.0))
          mesh.BoundList[iMarker][iElemSurface].changeOrientation();
      }
    }
  }
}

template <int DIM>
void FVMBuilder::SetEdges(MeshData<DIM> &mesh) {
  auto &pointList = mesh.PointList;
  Index nPoints = pointList.size();
  Index nEdge = 0, testEdge = 0;
  Index iPoint, jPoint, iNode, jNode;

  for (iPoint = 0; iPoint < nPoints; ++iPoint) {
    for (iNode = 0; iNode < pointList[iPoint].getnPoint(); ++iNode) {
      jPoint = pointList[iPoint].getPoint(iNode);
      for (jNode = 0; jNode < pointList[jPoint].getnPoint(); ++jNode) {
        if (pointList[jPoint].getPoint(jNode) == iPoint) {
          testEdge = pointList[jPoint].getEdge(jNode);
          break;
        }
      }
      if (testEdge == INVALID_INDEX) {
        pointList[iPoint].setEdge(iNode, nEdge);
        pointList[jPoint].setEdge(jNode, nEdge);
        nEdge++;
      }
    }
  }

  mesh.EdgeList.reserve(nEdge);

  for (iPoint = 0; iPoint < nPoints; ++iPoint) {
    for (iNode = 0; iNode < pointList[iPoint].getnPoint(); ++iNode) {
      jPoint = pointList[iPoint].getPoint(iNode);
      Index iEdge = findEdge(mesh, iPoint, jPoint);
      if (iPoint < jPoint) {
        mesh.EdgeList.emplace(mesh.EdgeList.begin() + iEdge, iPoint, jPoint);
      }
    }
  }
}

template <int DIM>
Index FVMBuilder::findEdge(const MeshData<DIM> &mesh, Index point1,
                           Index point2) {
  Index iPoint = 0;
  Index iNode = 0;
  auto &pointList = mesh.PointList;
  for (iNode = 0; iNode < pointList[point1].getnPoint(); ++iNode) {
    iPoint = pointList[point1].getPoint(iNode);
    if (iPoint == point2) {
      break;
    }
  }

  if (iPoint == point2) {
    return pointList[point1].getEdge(iNode);
  } else {
    std::cout << "\n\n   !!! Error !!!\n" << std::endl;
    std::cout << "Can't find the edge that connects " << point1 << " and "
              << point2 << "." << std::endl;
    return INVALID_INDEX;
  }
}

template <int DIM>
void FVMBuilder::SetVertex(MeshData<DIM> &mesh) {
  auto &pointList = mesh.PointList;
  Index nPoints = pointList.size();

  for (Index iPoint = 0; iPoint < nPoints; ++iPoint) {
    for (Index iMarker = 0; iMarker < mesh.nBound; ++iMarker) {
      pointList[iPoint].setVertex(iMarker, INVALID_INDEX);
    }
  }

  std::vector<Index> nVertex(mesh.nBound, 0);

  for (Index iMarker = 0; iMarker < mesh.nBound; ++iMarker) {
    for (Index iElem = 0; iElem < mesh.nElemBound[iMarker]; ++iElem) {
      for (Index iNode = 0; iNode < mesh.BoundList[iMarker][iElem].getnNodes();
           ++iNode) {
        Index iPoint = mesh.BoundList[iMarker][iElem].getNodes(iNode);

        if (pointList[iPoint].getVertex(iMarker) == INVALID_INDEX) {
          Index iVertex = nVertex[iMarker];
          pointList[iPoint].setVertex(iMarker, iVertex);
          nVertex[iMarker]++;
        }
      }
    }
  }

  for (Index iPoint = 0; iPoint < nPoints; ++iPoint) {
    for (Index iMarker = 0; iMarker < mesh.nBound; ++iMarker) {
      pointList[iPoint].setVertex(iMarker, INVALID_INDEX);
    }
  }

  mesh.VertexList.resize(mesh.nBound);

  for (Index iMarker = 0; iMarker < mesh.nBound; ++iMarker) {
    mesh.VertexList[iMarker].reserve(nVertex[iMarker]);
    nVertex[iMarker] = 0;

    for (Index iElem = 0; iElem < mesh.nElemBound[iMarker]; ++iElem) {
      for (Index iNode = 0; iNode < mesh.BoundList[iMarker][iElem].getnNodes();
           ++iNode) {
        Index iPoint = mesh.BoundList[iMarker][iElem].getNodes(iNode);

        if (pointList[iPoint].getVertex(iMarker) == INVALID_INDEX) {
          Index iVertex = nVertex[iMarker];
          mesh.VertexList[iMarker].push_back(iPoint);
          pointList[iPoint].setVertex(iMarker, iVertex);
          nVertex[iMarker]++;
        }
      }
    }
  }
}

template <int DIM>
void FVMBuilder::SetCG(MeshData<DIM> &mesh) {
  auto &pointList = mesh.PointList;
  auto &elemList = mesh.ElemList;
  auto &boundList = mesh.BoundList;
  auto &edgeList = mesh.EdgeList;

  Index nPoints = pointList.size();
  Index nElem = elemList.size();

  std::vector<std::array<real, DIM>> coord;

  for (Index iElem = 0; iElem < nElem; ++iElem) {
    Index nNode = elemList[iElem].getnNodes();
    coord.resize(nNode);

    for (Index iNode = 0; iNode < nNode; ++iNode) {
      Index elemPoint = elemList[iElem].getNodes(iNode);
      for (Index iDim = 0; iDim < DIM; ++iDim) {
        coord[iNode][iDim] = pointList[elemPoint].getCoord(iDim);
      }
    }

    for (Index iNode = 0; iNode < nNode; ++iNode) {
      for (Index iDim = 0; iDim < DIM; ++iDim) {
        elemList[iElem].setCoordCenter(iDim, coord[iNode][iDim] / real(nNode));
      }
    }

    for (Index iFace = 0; iFace < elemList[iElem].getnFaces(); ++iFace) {
      Index nNodesFace = elemList[iElem].getnNodesFace(iFace);
      for (Index iDim = 0; iDim < DIM; ++iDim) {
        for (Index iNode = 0; iNode < nNodesFace; ++iNode) {
          Index nodeFace = elemList[iElem].getFaces(iFace, iNode);
          elemList[iElem].setCoordFace(
              iFace, iDim, coord[nodeFace][iDim] / real(nNodesFace));
        }
      }
    }
  }

  for (Index iMarker = 0; iMarker < mesh.nBound; ++iMarker) {
    for (Index iElem = 0; iElem < mesh.nElemBound[iMarker]; ++iElem) {
      Index nNode = boundList[iMarker][iElem].getnNodes();
      coord.resize(nNode);

      for (Index iNode = 0; iNode < nNode; ++iNode) {
        Index elemPoint = boundList[iMarker][iElem].getNodes(iNode);
        for (Index iDim = 0; iDim < DIM; ++iDim) {
          coord[iNode][iDim] = pointList[elemPoint].getCoord(iDim);
        }
      }

      for (Index iNode = 0; iNode < nNode; ++iNode) {
        for (Index iDim = 0; iDim < DIM; ++iDim) {
          boundList[iMarker][iElem].setCoordCenter(
              iDim, coord[iNode][iDim] / real(nNode));
        }
      }

      for (Index iFace = 0; iFace < boundList[iMarker][iElem].getnFaces();
           ++iFace) {
        Index nNodesFace = boundList[iMarker][iElem].getnNodesFace(iFace);
        for (Index iDim = 0; iDim < DIM; ++iDim) {
          for (Index iNode = 0; iNode < nNodesFace; ++iNode) {
            Index nodeFace = boundList[iMarker][iElem].getFaces(iFace, iNode);
            boundList[iMarker][iElem].setCoordFace(
                iFace, iDim, coord[nodeFace][iDim] / real(nNodesFace));
          }
        }
      }
    }
  }

  std::array<std::array<real, DIM>, 2> edgeCoord{};
  for (Index iEdge = 0; iEdge < edgeList.size(); ++iEdge) {
    Index nNode = edgeList[iEdge].getnNodes();

    for (Index iNode = 0; iNode < nNode; ++iNode) {
      Index edgePoint = edgeList[iEdge].getNode(iNode);

      for (Index iDim = 0; iDim < DIM; ++iDim) {
        edgeCoord[iNode][iDim] = pointList[edgePoint].getCoord(iDim);
      }
    }

    edgeList[iEdge].setCoord(edgeCoord);
  }
}

template <int DIM>
void FVMBuilder::SetControlVolume(MeshData<DIM> &mesh) {
  auto &elemList = mesh.ElemList;
  auto &pointList = mesh.PointList;
  auto &edgeList = mesh.EdgeList;
  Index nPoints = pointList.size();
  Index nElem = elemList.size();
  Index nEdgesFace, faceiPoint, facejPoint;
  real domainVolume = 0.0;
  real TotalDomainVolume = 0.0;

  std::array<real, DIM> coordEdge, coordFaceElem, coordElem, coordFaceiPoint,
      coordFacejPoint, normalFace;

  for (Index iElem = 0; iElem < nElem; ++iElem) {
    for (Index iFace = 0; iFace < elemList[iElem].getnFaces(); ++iFace) {
      if (DIM == 2) {
        nEdgesFace = 1;
      }
      if (DIM == 3) {
        nEdgesFace = elemList[iElem].getnNodesFace(iFace);
      }

      for (Index iEdgesFace = 0; iEdgesFace < nEdgesFace; ++iEdgesFace) {
        if (DIM == 2) {
          faceiPoint =
              elemList[iElem].getNodes(elemList[iElem].getFaces(iFace, 0));
          facejPoint =
              elemList[iElem].getNodes(elemList[iElem].getFaces(iFace, 1));
        }

        if (DIM == 3) {
          faceiPoint = elemList[iElem].getNodes(
              elemList[iElem].getFaces(iFace, iEdgesFace));
          if (iEdgesFace != nEdgesFace - 1) {
            facejPoint = elemList[iElem].getNodes(
                elemList[iElem].getFaces(iFace, iEdgesFace + 1));
          } else {
            facejPoint =
                elemList[iElem].getNodes(elemList[iElem].getFaces(iFace, 0));
          }
        }

        bool changeFaceOrientation = false;
        if (faceiPoint > facejPoint) {
          changeFaceOrientation = true;
        }
        Index iEdge = findEdge(mesh, faceiPoint, facejPoint);

        for (Index iDim = 0; iDim < DIM; ++iDim) {
          coordEdge[iDim] = edgeList[iEdge].getCoord(iDim);
          coordElem[iDim] = elemList[iElem].getCoordCenter(iDim);
          coordFaceElem[iDim] = elemList[iElem].getCoordFace(iFace, iDim);
          coordFaceiPoint[iDim] = pointList[faceiPoint].getCoord(iDim);
          coordFacejPoint[iDim] = pointList[facejPoint].getCoord(iDim);
        }

        if (DIM == 2) {
          if (changeFaceOrientation) {
            edgeList[iEdge].setNodesCoord(coordElem, coordEdge);
          } else {
            edgeList[iEdge].setNodesCoord(coordEdge, coordElem);
          }
          real area =
              edgeList[iEdge].getVolume(coordFaceiPoint, coordEdge, coordElem);
          pointList[faceiPoint].addVolume(area);
          domainVolume += area;
          area =
              edgeList[iEdge].getVolume(coordFacejPoint, coordEdge, coordElem);
          pointList[facejPoint].addVolume(area);
          domainVolume += area;
        } else if (DIM == 3) {
          if (changeFaceOrientation) {
            edgeList[iEdge].setNodesCoord(coordFaceElem, coordEdge, coordElem);
          } else {
            edgeList[iEdge].setNodesCoord(coordEdge, coordFaceElem, coordElem);
          }
          real volume = edgeList[iEdge].getVolume(coordFaceiPoint, coordEdge,
                                                  coordFaceElem, coordElem);
          pointList[faceiPoint].addVolume(volume);
          domainVolume += volume;
          volume = edgeList[iEdge].getVolume(coordFacejPoint, coordEdge,
                                             coordFaceElem, coordElem);
          pointList[facejPoint].addVolume(volume);
          domainVolume += volume;
        }
      }
    }
  }

  for (Index iEdge = 0; iEdge < edgeList.size(); ++iEdge) {
    normalFace = edgeList[iEdge].getNormal();
    real area = 0.0;
    for (Index iDim = 0; iDim < DIM; ++iDim) {
      area += normalFace[iDim] * normalFace[iDim];
    }
    area = std::sqrt(area);

    if (area == 0.0) {
      for (Index iDim = 0; iDim < DIM; ++iDim) {
        normalFace[iDim] = std::numeric_limits<real>::epsilon() *
                           std::numeric_limits<real>::epsilon();
      }
      edgeList[iEdge].setNormal(normalFace);
    }
  }

  TotalDomainVolume = domainVolume;

  if (DIM == 2)
    std::cout << "Area of the computational grid: " << TotalDomainVolume << "."
              << std::endl;
  if (DIM == 3)
    std::cout << "Volume of the computational grid: " << TotalDomainVolume
              << "." << std::endl;
}

template <int DIM>
void FVMBuilder::SetBoundControlVolume(MeshData<DIM> &mesh) {
  auto &edgeList = mesh.EdgeList;
  auto &pointList = mesh.PointList;
  auto &boundList = mesh.BoundList;
  auto &vertexList = mesh.VertexList;
  Index nMarker = mesh.nBound;
  std::array<real, DIM> coordEdge, coordElem, coordVertex;

  for (Index iMarker = 0; iMarker < nMarker; ++iMarker) {
    for (Index iElem = 0; iElem < mesh.nElemBound[iMarker]; ++iElem) {
      for (Index iNode = 0; iNode < boundList[iMarker][iElem].getnNodes();
           ++iNode) {
        Index iPoint = boundList[iMarker][iElem].getNodes(iNode);
        Index iVertex = pointList[iPoint].getVertex(iMarker);

        for (Index iNeighborNodes = 0;
             iNeighborNodes <
             boundList[iMarker][iElem].getnNeighborNodes(iNode);
             ++iNeighborNodes) {
          Index neighborNode =
              boundList[iMarker][iElem].getNodesNeighbor(iNode, iNeighborNodes);
          Index neighborPoint =
              boundList[iMarker][iElem].getNodes(neighborNode);

          Index iEdge = findEdge(mesh, iPoint, neighborPoint);

          for (Index iDim = 0; iDim < DIM; ++iDim) {
            coordEdge[iDim] = edgeList[iEdge].getCoord(iDim);
            coordElem[iDim] = boundList[iMarker][iElem].getCoordCenter(iDim);
            coordVertex[iDim] = pointList[iPoint].getCoord(iDim);
          }

          if (DIM == 2) {
            if (iNode == 0) {
              vertexList[iMarker][iVertex].setNodesCoord(coordElem,
                                                         coordVertex);
            }
            if (iNode == 1) {
              vertexList[iMarker][iVertex].setNodesCoord(coordVertex,
                                                         coordElem);
            }
          } else if (DIM == 3) {
            if (iNeighborNodes == 0) {
              vertexList[iMarker][iVertex].setNodesCoord(coordElem, coordEdge,
                                                         coordVertex);
            }
            if (iNeighborNodes == 1) {
              vertexList[iMarker][iVertex].setNodesCoord(coordEdge, coordElem,
                                                         coordVertex);
            }
          }
        }
      }
    }
  }

  for (Index iMarker = 0; iMarker < nMarker; ++iMarker) {
    for (Index iVertex = 0; iVertex < mesh.VertexList[iMarker].size();
         ++iVertex) {
      std::array<real, DIM> normalFace =
          vertexList[iMarker][iVertex].getNormal();
      real area = 0.0;
      for (Index iDim = 0; iDim < DIM; ++iDim) {
        area += normalFace[iDim] * normalFace[iDim];
      }
      area = std::sqrt(area);
      if (area == 0.0) {
        for (Index iDim = 0; iDim < DIM; ++iDim) {
          normalFace[iDim] = std::numeric_limits<real>::epsilon() *
                             std::numeric_limits<real>::epsilon();
        }
        vertexList[iMarker][iVertex].setNormal(normalFace);
      }
    }
  }
}

template <int DIM>
void FVMBuilder::FindNormalNeighbor(MeshData<DIM> &mesh) {
  auto &edgeList = mesh.EdgeList;
  auto &pointList = mesh.PointList;
  auto &boundList = mesh.BoundList;
  auto &vertexList = mesh.VertexList;
  Index nMarker = mesh.nBound;

  for (Index iMarker = 0; iMarker < nMarker; ++iMarker) {
    for (Index iVertex = 0; iVertex < vertexList[iMarker].size(); ++iVertex) {
      Index iPoint = vertexList[iMarker][iVertex].getNode();
      std::array<real, DIM> normal = vertexList[iMarker][iVertex].getNormal();

      real pointNormal = 0.0;
      real cosMax = -1.0;
      for (Index iNeighbor = 0; iNeighbor < pointList[iPoint].getnPoint();
           ++iNeighbor) {
        Index jPoint = pointList[iPoint].getPoint(iNeighbor);
        real scalarProd = 0.0;
        real normVect = 0.0;
        real normNormal = 0.0;
        for (Index iDim = 0; iDim < DIM; ++iDim) {
          real diffCoord = pointList[jPoint].getCoord(iDim) -
                           pointList[iPoint].getCoord(iDim);
          scalarProd += diffCoord * normal[iDim];
          normVect += diffCoord * diffCoord;
          normNormal += normal[iDim] * normal[iDim];
        }

        normVect = std::sqrt(normVect);
        normNormal = std::sqrt(normNormal);
        real cosAlpha = scalarProd / (normVect * normNormal);

        if (cosAlpha >= cosMax) {
          pointNormal = jPoint;
          cosMax = cosAlpha;
        }
      }

      vertexList[iMarker][iVertex].setNormalNeighbor(pointNormal);
    }
  }
}
}  // namespace GridBuilder
}  // namespace preprocess