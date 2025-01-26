#pragma once
#include <pre/macro.h>
#include <pre/reader.h>

#include <string>
namespace preprocess {

class Geometry {
public:
  int totNodes, phyNodes, totEdges, phyEdges, numTria, numBoundSegs,
      numBoundFaces, numBoundNodes;

  int *BoundTypes;

  Tria *tria;

  Edge *edge;

  BoundaryFace *boundaryFace;

  BoundaryNode *boundaryNode;

  idBoundary *ibound;

  Node *coords;
  Node *sij;

  double *vol;

  Node *sbf;
  Node *sproj;

  std::string *bname;

  Geometry(const std::string &filename, const std::string &commentChar = "#");
  ~Geometry();
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
  EdgeList *tmpElist;

  void CheckMetrics();
  void DeleteTmpElist();
  void DummyNodes();
  void FaceVectorsSymm();
  void FaceVectorsVolumes();
  void FaceVectorsVolumesBound();
  void volumeProjections();
};
} // namespace preprocess