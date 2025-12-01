#pragma once

#include <element/element.hpp>
#include <grid/grid.hpp>
#include <mesh/mesh_data.hpp>

namespace preprocess {
namespace GridBuilder {

struct FVMBuilder {
  static void build(MeshDataBase &mesh);

  template <int DIM>
  static void buildImpl(MeshData<DIM> &mesh);

  template <int DIM>
  static void SetPointConnectivity(MeshData<DIM> &mesh);

  template <int DIM>
  static void SetRCMOrdering(MeshData<DIM> &mesh);

  template <int DIM>
  static void SetElementConnectivity(MeshData<DIM> &mesh);

  template <int DIM>
  static bool findFace(MeshData<DIM> &mesh, Index elem1, Index elem2,
                       Index &elem1Face, Index &elem2Face);

  template <int DIM>
  static void SetBoundVolume(MeshData<DIM> &mesh);

  template <int DIM>
  static void CheckOrientation(MeshData<DIM> &mesh);

  template <int DIM>
  static void SetEdges(MeshData<DIM> &mesh);

  template <int DIM>
  static Index findEdge(const MeshData<DIM> &mesh, Index point1, Index point2);

  template <int DIM>
  static void SetVertex(MeshData<DIM> &mesh);

  template <int DIM>
  static void SetCG(MeshData<DIM> &mesh);

  template <int DIM>
  static void SetControlVolume(MeshData<DIM> &mesh);

  template <int DIM>
  static void SetBoundControlVolume(MeshData<DIM> &mesh);

  template <int DIM>
  static void FindNormalNeighbor(MeshData<DIM> &mesh);
};
}  // namespace GridBuilder
}  // namespace preprocess
