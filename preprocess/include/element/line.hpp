#pragma once

#include <cassert>
#include <element/element_base.hpp>
namespace preprocess {
template <int DIM>
class Line : public ElementBase<DIM, Line<DIM>, LineNodes, LineFaces> {
 public:
  using BaseType = ElementBase<DIM, Line<DIM>, LineNodes, LineFaces>;
  Line() = default;
  Line(Index point_0, Index point_1);
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
  constexpr static short nFaces = 1;
  constexpr static short nNodes = 2;
  constexpr static short nNeighbors = 1;
  constexpr static short VTKtype = 3;
  constexpr static short maxNodesFace = 2;

  constexpr static short Faces[1][2] = {{0, 1}};
  constexpr static short NodesNeighbor[2][1] = {{1}, {0}};
  constexpr static short nNodesFace[1] = {2};
  constexpr static short nNeighborNodes[2] = {1, 1};

  std::array<real, DIM> coordCenter{};
  std::array<std::array<real, DIM>, nFaces> coordFaces{};
};

template <int DIM>
constexpr Index Line<DIM>::getMaxNodesFaceImpl() noexcept {
  return maxNodesFace;
}

template <int DIM>
constexpr Index Line<DIM>::getVTKtypeImpl() noexcept {
  return VTKtype;
}

template <int DIM>
Line<DIM>::Line(Index point_0, Index point_1) {
  this->nodes[0] = point_0;
  this->nodes[1] = point_1;
  for (Index iDim = 0; iDim < DIM; ++iDim) {
    coordCenter[iDim] = 0.0;
    for (Index iface = 0; iface < nFaces; ++iface) {
      coordFaces[iface][iDim] = 0.0;
    }
  }

  this->neighborElements.fill(INVALID_INDEX);
}

template <int DIM>
real Line<DIM>::getCoordCenterImpl(Index idim) const noexcept {
  assert(idim >= 0 && idim < DIM && "index out of range");
  return coordCenter[idim];
}

template <int DIM>
void Line<DIM>::setCoordCenterImpl(Index idim, real icoord) noexcept {
  assert(idim >= 0 && idim < DIM && "index out of range");
  coordCenter[idim] += icoord;
}

template <int DIM>
real Line<DIM>::getCoordFaceImpl(Index iface, Index idim) const noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(idim >= 0 && idim < DIM && "index out of range");
  return coordFaces[iface][idim];
}

template <int DIM>
void Line<DIM>::setCoordFaceImpl(Index iface, Index idim,
                                 real icoord) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(idim >= 0 && idim < DIM && "index out of range");
  coordFaces[iface][idim] += icoord;
}

template <int DIM>
Index Line<DIM>::getNodesNeighborImpl(Index inode, Index ineighbor) noexcept {
  assert(inode >= 0 && inode < nNodes && "index out of range");
  assert(ineighbor >= 0 && ineighbor < nNeighborNodes[inode] &&
         "index out of range");
  return NodesNeighbor[inode][ineighbor];
}

template <int DIM>
Index Line<DIM>::getnNodesFaceImpl(Index iface) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  return nNodesFace[iface];
}

template <int DIM>
Index Line<DIM>::getnNeighborNodesImpl(Index inode) noexcept {
  assert(inode >= 0 && inode < nNodes && "index out of range");
  return nNeighborNodes[inode];
}

template <int DIM>
Index Line<DIM>::getFacesImpl(Index iface, Index inode) noexcept {
  assert(iface >= 0 && iface < nFaces && "index out of range");
  assert(inode >= 0 && inode < nNodesFace[iface] && "index out of range");
  return Faces[iface][inode];
}

template <int DIM>
void Line<DIM>::changeOrientationImpl() noexcept {
  std::swap(this->nodes[0], this->nodes[1]);
}
}  // namespace preprocess