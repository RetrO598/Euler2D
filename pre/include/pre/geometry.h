#pragma once
#include <pre/macro.h>
#include <pre/reader.h>

#include <string>
#include <vector>
namespace preprocess {

class Geometry {
 public:
  int phyNodes, phyEdges, numTria, numBoundSegs, numBoundFaces, numBoundNodes;

  std::vector<int> BoundTypes;
  std::vector<Tria> tria;
  std::vector<Edge> edge;
  std::vector<BoundaryFace> boundaryFace;
  std::vector<BoundaryNode> boundaryNode;
  std::vector<idBoundary> ibound;
  std::vector<Node> coords;
  std::vector<Node> sij;
  std::vector<double> vol;
  std::vector<Node> sbf;
  std::vector<Node> sproj;
  std::vector<std::string> bname;
  std::vector<vertex> vertexList;

  Geometry(const std::string &filename, const std::string &commentChar = "#");
  Geometry(const Geometry &geomtry) = delete;
  Geometry(Geometry &&geometry) = delete;
  Geometry &operator=(const Geometry &geomtry) = delete;
  Geometry &&operator=(Geometry &&geometry) = delete;

  void ComputeMetrics();
  void GenerateEdgeList();
  void printInfo();
  void GetNumberBoundNodes(int btypeFrom, int btypeTo) const;
  void ReadGrid();
  void outputMeshInfo();

 private:
  simpleReader gridReader;
  std::vector<EdgeList> tmpElist;

  void CheckMetrics();
  // void DeleteTmpElist();
  void DummyNodes();
  void FaceVectorsSymm();
  void FaceVectorsVolumes();
  void FaceVectorsVolumesBound();
  void volumeProjections();
};
}  // namespace preprocess