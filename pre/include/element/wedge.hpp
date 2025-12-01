#pragma once
#include <cassert>
#include <element/element_base.hpp>
#include <element/utils.hpp>
namespace preprocess {
template <int DIM>
class Wedge : public ElementBase<DIM, Wedge<DIM>, WedgeNodes, WedgeFaces> {
 public:
  using BaseType = ElementBase<DIM, Wedge<DIM>, WedgeNodes, WedgeFaces>;
  Wedge() = default;
  Wedge(Index point_0, Index point_1, Index point_2, Index point_3,
        Index point_4, Index point_5);
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
  constexpr static short nNodes = 6;
  constexpr static short nNeighbors = 5;
  constexpr static short VTKtype = 13;
  constexpr static short maxNodesFace = 4;

  constexpr static short Faces[5][4] = {
      {3, 4, 1, 0}, {5, 2, 1, 4}, {2, 5, 3, 0}, {0, 1, 2, 2}, {5, 4, 3, 3}};

  constexpr static short NodesNeighbor[6][3] = {
      {1, 2, 3}, {0, 2, 4}, {1, 0, 5}, {0, 4, 5}, {3, 5, 1}, {4, 3, 2}};

  constexpr static short nNodesFace[5] = {4, 4, 4, 3, 3};
  constexpr static short nNeighborNodes[6] = {3, 3, 3, 3, 3, 3};

  std::array<real, DIM> coordCenter{};
  std::array<std::array<real, DIM>, nFaces> coordFaces{};
};

template <int DIM>
constexpr Index Wedge<DIM>::getMaxNodesFaceImpl() noexcept {
  return maxNodesFace;
}

template <int DIM>
constexpr Index Wedge<DIM>::getVTKtypeImpl() noexcept {
  return VTKtype;
}

template <int DIM>
Wedge<DIM>::Wedge(Index point_0, Index point_1, Index point_2, Index point_3,
                  Index point_4, Index point_5) {
  this->nodes[0] = point_0;
  this->nodes[1] = point_1;
  this->nodes[2] = point_2;
  this->nodes[3] = point_3;
  this->nodes[4] = point_4;
  this->nodes[5] = point_5;

  for (Index iDim = 0; iDim < DIM; ++iDim) {
    coordCenter[iDim] = 0.0;
    for (Index iface = 0; iface < nFaces; ++iface) {
      coordFaces[iface][iDim] = 0.0;
    }
  }

  this->neighborElements.fill(INVALID_INDEX);
}

template <int DIM>
real Wedge<DIM>::getCoordCenterImpl(Index idim) const noexcept {
  assert(idim >= 0 && idim < DIM && "index out of range");
  return coordCenter[idim];
}

template <int DIM>
void Wedge<DIM>::setCoordCenterImpl(Index idim, real icoord) noexcept {
  assert(idim >= 0 && idim < DIM && "index out of range");
  coordCenter[idim] += icoord;
}

template <int DIM>
real Wedge<DIM>::getCoordFaceImpl(Index iface, Index idim) const noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(idim >= 0 && idim < DIM && "index out of range");
  return coordFaces[iface][idim];
}

template <int DIM>
void Wedge<DIM>::setCoordFaceImpl(Index iface, Index idim,
                                  real icoord) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(idim >= 0 && idim < DIM && "index out of range");
  coordFaces[iface][idim] += icoord;
}

template <int DIM>
Index Wedge<DIM>::getNodesNeighborImpl(Index inode, Index ineighbor) noexcept {
  assert(inode >= 0 && inode < nNodes && "index out of range");
  assert(ineighbor >= 0 && ineighbor < nNeighborNodes[inode] &&
         "index out of range");
  return NodesNeighbor[inode][ineighbor];
}

template <int DIM>
Index Wedge<DIM>::getnNodesFaceImpl(Index iface) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  return nNodesFace[iface];
}

template <int DIM>
Index Wedge<DIM>::getnNeighborNodesImpl(Index inode) noexcept {
  assert(inode >= 0 && inode < nNodes && "index out of range");
  return nNeighborNodes[inode];
}

template <int DIM>
Index Wedge<DIM>::getFacesImpl(Index iface, Index inode) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(inode >= 0 && inode < nNodesFace[iface] && "index out of range");
  return Faces[iface][inode];
}

template <int DIM>
void Wedge<DIM>::changeOrientationImpl() noexcept {
  std::swap(this->nodes[0], this->nodes[1]);
  std::swap(this->nodes[3], this->nodes[4]);
}

}  // namespace preprocess