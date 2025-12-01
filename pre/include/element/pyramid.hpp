#pragma once
#include <cassert>
#include <element/element_base.hpp>
#include <element/utils.hpp>
namespace preprocess {
template <int DIM>
class Pyramid
    : public ElementBase<DIM, Pyramid<DIM>, PyramidNodes, PyramidFaces> {
 public:
  using BaseType = ElementBase<DIM, Pyramid<DIM>, PyramidNodes, PyramidFaces>;
  Pyramid() = default;
  Pyramid(Index point_0, Index point_1, Index point_2, Index point_3,
          Index point_4);
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
  constexpr static short nFaces = 5;
  constexpr static short nNodes = 5;
  constexpr static short nNeighbors = 5;
  constexpr static short VTKtype = 14;
  constexpr static short maxNodesFace = 4;

  constexpr static short Faces[5][4] = {
      {0, 3, 2, 1}, {4, 3, 0, 0}, {4, 0, 1, 1}, {2, 4, 1, 1}, {3, 4, 2, 2}};

  constexpr static short NodesNeighbor[5][4] = {
      {1, 3, 4, 4}, {0, 2, 4, 4}, {1, 3, 4, 4}, {2, 0, 4, 4}, {0, 1, 2, 3}};

  constexpr static short nNodesFace[5] = {4, 3, 3, 3, 3};

  constexpr static short nNeighborNodes[5] = {3, 3, 3, 3, 4};

  std::array<real, DIM> coordCenter{};
  std::array<std::array<real, DIM>, nFaces> coordFaces{};
};

template <int DIM>
constexpr Index Pyramid<DIM>::getMaxNodesFaceImpl() noexcept {
  return maxNodesFace;
}

template <int DIM>
constexpr Index Pyramid<DIM>::getVTKtypeImpl() noexcept {
  return VTKtype;
}

template <int DIM>
Pyramid<DIM>::Pyramid(Index point_0, Index point_1, Index point_2,
                      Index point_3, Index point_4) {
  this->nodes[0] = point_0;
  this->nodes[1] = point_1;
  this->nodes[2] = point_2;
  this->nodes[3] = point_3;
  this->nodes[4] = point_4;
  for (Index iDim = 0; iDim < DIM; ++iDim) {
    coordCenter[iDim] = 0.0;
    for (Index iface = 0; iface < nFaces; ++iface) {
      coordFaces[iface][iDim] = 0.0;
    }
  }
  this->neighborElements.fill(INVALID_INDEX);
}

template <int DIM>
real Pyramid<DIM>::getCoordCenterImpl(Index idim) const noexcept {
  assert(idim >= 0 && idim < DIM && "index out of range");
  return coordCenter[idim];
}

template <int DIM>
void Pyramid<DIM>::setCoordCenterImpl(Index idim, real icoord) noexcept {
  assert(idim >= 0 && idim < DIM && "index out of range");
  coordCenter[idim] += icoord;
}

template <int DIM>
real Pyramid<DIM>::getCoordFaceImpl(Index iface, Index idim) const noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(idim >= 0 && idim < DIM && "index out of range");
  return coordFaces[iface][idim];
}

template <int DIM>
void Pyramid<DIM>::setCoordFaceImpl(Index iface, Index idim,
                                    real icoord) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(idim >= 0 && idim < DIM && "index out of range");
  coordFaces[iface][idim] += icoord;
}

template <int DIM>
Index Pyramid<DIM>::getNodesNeighborImpl(Index inode,
                                         Index ineighbor) noexcept {
  assert(inode >= 0 && inode < nNodes && "index out of range");
  assert(ineighbor >= 0 && ineighbor < nNeighborNodes[inode] &&
         "index out of range");
  return NodesNeighbor[inode][ineighbor];
}

template <int DIM>
Index Pyramid<DIM>::getnNodesFaceImpl(Index iface) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  return nNodesFace[iface];
}

template <int DIM>
Index Pyramid<DIM>::getnNeighborNodesImpl(Index inode) noexcept {
  assert(inode >= 0 && inode < nNodes && "index out of range");
  return nNeighborNodes[inode];
}

template <int DIM>
Index Pyramid<DIM>::getFacesImpl(Index iface, Index inode) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(inode >= 0 && inode < nNodesFace[iface] && "index out of range");
  return Faces[iface][inode];
}

template <int DIM>
void Pyramid<DIM>::changeOrientationImpl() noexcept {
  std::swap(this->nodes[1], this->nodes[3]);
}
}  // namespace preprocess