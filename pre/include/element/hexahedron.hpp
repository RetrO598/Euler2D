#pragma once
#include <array>
#include <cassert>
#include <element/element_base.hpp>
#include <element/utils.hpp>

#include "common/config.hpp"

namespace preprocess {
template <int DIM>
class Hexahedron : public ElementBase<DIM, Hexahedron<DIM>, HexahedronNodes,
                                      HexahedronFaces> {
 public:
  using BaseType =
      ElementBase<DIM, Hexahedron<DIM>, HexahedronNodes, HexahedronFaces>;
  Hexahedron() = default;
  Hexahedron(Index point_0, Index point_1, Index point_2, Index point_3,
             Index point_4, Index point_5, Index point_6, Index point_7);

  constexpr static Index getMaxNodesFaceImpl() noexcept;

  constexpr static Index getVTKtypeImpl() noexcept;

  real getCoordCenterImpl(Index idim) const noexcept;

  void setCoordCenterImpl(Index idim, real icoord) noexcept;

  real getCoordFaceImpl(Index iface, Index idim) const noexcept;

  void setCoordFaceImpl(Index iface, Index idim, real icoord) noexcept;

  Index static getNodesNeighborImpl(Index inode, Index ineighbor) noexcept;

  Index static getnNodesFaceImpl(Index iface) noexcept;

  Index static getnNeighborNodesImpl(Index inode) noexcept;

  Index static getFacesImpl(Index iface, Index inode) noexcept;

  void changeOrientationImpl() noexcept;

 private:
  constexpr static short nFaces = 6;
  constexpr static short nNodes = 8;
  constexpr static short nNeighbors = 6;
  constexpr static short VTKtype = 12;
  constexpr static short maxNodesFace = 4;

  constexpr static short Faces[6][4] = {{0, 1, 5, 4}, {1, 2, 6, 5},
                                        {2, 3, 7, 6}, {3, 0, 4, 7},
                                        {0, 3, 2, 1}, {4, 5, 6, 7}};

  constexpr static short NodesNeighbor[8][3] = {{1, 3, 4}, {0, 2, 5}, {1, 3, 6},
                                                {0, 2, 7}, {0, 5, 7}, {4, 6, 1},
                                                {2, 5, 7}, {4, 3, 6}};

  constexpr static short nNodesFace[6] = {4, 4, 4, 4, 4, 4};
  constexpr static short nNeighborNodes[8] = {3, 3, 3, 3, 3, 3, 3, 3};

  std::array<real, DIM> coordCenter{};
  std::array<std::array<real, DIM>, nFaces> coordFaces{};
};

template <int DIM>
constexpr Index Hexahedron<DIM>::getMaxNodesFaceImpl() noexcept {
  return maxNodesFace;
}

template <int DIM>
constexpr Index Hexahedron<DIM>::getVTKtypeImpl() noexcept {
  return VTKtype;
}

template <int DIM>
Hexahedron<DIM>::Hexahedron(Index point_0, Index point_1, Index point_2,
                            Index point_3, Index point_4, Index point_5,
                            Index point_6, Index point_7) {
  this->nodes[0] = point_0;
  this->nodes[1] = point_1;
  this->nodes[2] = point_2;
  this->nodes[3] = point_3;
  this->nodes[4] = point_4;
  this->nodes[5] = point_5;
  this->nodes[6] = point_6;
  this->nodes[7] = point_7;
  for (Index iDim = 0; iDim < DIM; ++iDim) {
    coordCenter[iDim] = 0.0;
    for (Index iface = 0; iface < nFaces; ++iface) {
      coordFaces[iface][iDim] = 0.0;
    }
  }

  this->neighborElements.fill(INVALID_INDEX);
}

template <int DIM>
real Hexahedron<DIM>::getCoordCenterImpl(Index idim) const noexcept {
  assert(idim >= 0 && idim < DIM && "index out of range");
  return coordCenter[idim];
}

template <int DIM>
void Hexahedron<DIM>::setCoordCenterImpl(Index idim, real icoord) noexcept {
  assert(idim >= 0 && idim < DIM && "index out of range");
  coordCenter[idim] += icoord;
}

template <int DIM>
real Hexahedron<DIM>::getCoordFaceImpl(Index iface, Index idim) const noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(idim >= 0 && idim < DIM && "index out of range");
  return coordFaces[iface][idim];
}

template <int DIM>
void Hexahedron<DIM>::setCoordFaceImpl(Index iface, Index idim,
                                       real icoord) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(idim >= 0 && idim < DIM && "index out of range");
  coordFaces[iface][idim] += icoord;
}

template <int DIM>
Index Hexahedron<DIM>::getNodesNeighborImpl(Index inode,
                                            Index ineighbor) noexcept {
  assert(inode >= 0 && inode < nNodes && "index out of range");
  assert(ineighbor >= 0 && ineighbor < nNeighborNodes[inode] &&
         "index out of range");
  return NodesNeighbor[inode][ineighbor];
}

template <int DIM>
Index Hexahedron<DIM>::getnNodesFaceImpl(Index iface) noexcept {
  assert(iface >= 0 && iface < nFaces < DIM && "index out of range");
  return nNodesFace[iface];
}

template <int DIM>
Index Hexahedron<DIM>::getnNeighborNodesImpl(Index inode) noexcept {
  return nNeighborNodes[inode];
}

template <int DIM>
Index Hexahedron<DIM>::getFacesImpl(Index iface, Index inode) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(inode >= 0 && inode < nNodesFace[iface] && "index out of range");
  return Faces[iface][inode];
}

template <int DIM>
void Hexahedron<DIM>::changeOrientationImpl() noexcept {
  auto point0 = this->nodes[0];
  auto point1 = this->nodes[1];
  auto point2 = this->nodes[2];
  auto point3 = this->nodes[3];
  auto point4 = this->nodes[4];
  auto point5 = this->nodes[5];
  auto point6 = this->nodes[6];
  auto point7 = this->nodes[7];

  this->nodes[0] = point7;
  this->nodes[1] = point4;
  this->nodes[2] = point5;
  this->nodes[3] = point6;
  this->nodes[4] = point3;
  this->nodes[5] = point0;
  this->nodes[6] = point1;
  this->nodes[7] = point2;
}
}  // namespace preprocess