#pragma once
#include <pre/macro.h>
#include <pre/reader.h>

#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "element/element_concept.hpp"
#include "grid/point.hpp"
#include "pre/parameter.h"
namespace preprocess {

class Geometry {
 public:
  int phyNodes, phyEdges, numElems, numBoundSegs, numBoundFaces, numBoundNodes;

  std::vector<ElementHandle> elements;
  std::vector<Point<2>> pointList;
  std::vector<Edge2D> edge;
  std::vector<BoundaryFace> boundaryFace;
  std::vector<idBoundary> ibound;
  std::vector<Node> coords;
  std::vector<Node> sij;
  std::vector<double> vol;
  std::vector<Node> sbf;
  std::vector<Node> sproj;
  std::vector<std::string> bname;
  std::vector<vertex> vertexList;
  std::unordered_map<std::string, preprocess::BoundaryType> boundaryMap;

  Geometry() = default;
  Geometry(const parameter &param, const std::string &commentChar = "#");
  Geometry(const Geometry &geomtry) = delete;
  Geometry(Geometry &&geometry) {
    this->bname = std::move(geometry.bname);
    this->boundaryFace = std::move(geometry.boundaryFace);
    this->boundaryMap = std::move(geometry.boundaryMap);
    this->numBoundFaces = geometry.numBoundFaces;
    this->numBoundNodes = geometry.numBoundNodes;
    this->numBoundSegs = geometry.numBoundSegs;
    this->ibound = std::move(geometry.ibound);
    this->coords = std::move(geometry.coords);
    this->edge = std::move(geometry.edge);
    this->phyNodes = geometry.phyNodes;
    this->phyEdges = geometry.phyEdges;
    this->sbf = std::move(geometry.sbf);
    this->sij = std::move(geometry.sij);
    this->sproj = std::move(geometry.sproj);
    this->vertexList = std::move(geometry.vertexList);
    this->vol = std::move(geometry.vol);
    this->numElems = geometry.numElems;
    this->elements = std::move(geometry.elements);
    this->pointList = std::move(geometry.pointList);
    this->periodicInfo = std::move(geometry.periodicInfo);
    this->periodicMaster = std::move(geometry.periodicMaster);
  };
  Geometry &operator=(const Geometry &geomtry) = delete;
  Geometry &operator=(Geometry &&geometry) {
    this->bname = std::move(geometry.bname);
    this->boundaryFace = std::move(geometry.boundaryFace);
    this->boundaryMap = std::move(geometry.boundaryMap);
    this->numBoundFaces = geometry.numBoundFaces;
    this->numBoundNodes = geometry.numBoundNodes;
    this->numBoundSegs = geometry.numBoundSegs;
    this->ibound = std::move(geometry.ibound);
    this->coords = std::move(geometry.coords);
    this->edge = std::move(geometry.edge);
    this->phyNodes = geometry.phyNodes;
    this->phyEdges = geometry.phyEdges;
    this->sbf = std::move(geometry.sbf);
    this->sij = std::move(geometry.sij);
    this->sproj = std::move(geometry.sproj);
    this->vertexList = std::move(geometry.vertexList);
    this->vol = std::move(geometry.vol);
    this->numElems = geometry.numElems;
    this->elements = std::move(geometry.elements);
    this->pointList = std::move(geometry.pointList);
    this->periodicInfo = std::move(geometry.periodicInfo);
    this->periodicMaster = std::move(geometry.periodicMaster);

    return *this;
  };

  void ComputeMetrics();
  void printInfo();
  void outputMeshInfo();

  //  private:
  std::vector<EdgeList> tmpElist;
  std::vector<PeriodicInfo> periodicInfo;
  std::set<std::string> periodicMaster;

  void CheckMetrics();
  void FaceVectorsSymm();
  void volumeProjections();
  void ComputeWallDistance();
  void MatchPeriodic();
};
}  // namespace preprocess