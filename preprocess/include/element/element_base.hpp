#pragma once
#include <array>
#include <common/common.hpp>
#include <common/config.hpp>
#include <element/utils.hpp>

namespace preprocess {
template <int DIM, class Derived, Index Nodes, Index Faces>
class ElementBase {
 public:
  using ElementType = Derived;

  constexpr ElementType& crtpCast() noexcept;

  constexpr const ElementType& crtpCast() const noexcept;

  constexpr static Index getnNodes() noexcept;

  constexpr static Index getnFaces() noexcept;

  constexpr static Index getnNeighbors() noexcept;

  Index getNeighborElems(Index iface) const noexcept;

  void setNeighborElems(Index iface, Index elemIndex) noexcept;

  Index getNodes(Index inode) const noexcept;

  void setNodes(Index inode, Index nodeIndex) noexcept;

  constexpr Index getMaxNodesFace() const noexcept;

  constexpr Index getVTKtype() const noexcept;

  real getCoordCenter(Index idim) const noexcept;

  void setCoordCenter(Index idim, real icoord) noexcept;

  real getCoordFace(Index iface, Index idim) const noexcept;

  void setCoordFace(Index iface, Index idim, real icoord) noexcept;

  Index getNodesNeighbor(Index inode, Index ineighbor) const noexcept;

  Index getnNodesFace(Index iface) const noexcept;

  Index getnNeighborNodes(Index inode) const noexcept;

  void changeOrientation() noexcept;

  Index getFaces(Index iface, Index inode) const noexcept;

  std::array<Index, Nodes> nodes{};
  std::array<Index, Faces> neighborElements{};
};

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE constexpr Derived&
ElementBase<DIM, Derived, Nodes, Faces>::crtpCast() noexcept {
  return static_cast<Derived&>(*this);
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE constexpr const Derived&
ElementBase<DIM, Derived, Nodes, Faces>::crtpCast() const noexcept {
  return static_cast<const Derived&>(*this);
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE constexpr Index
ElementBase<DIM, Derived, Nodes, Faces>::getnNodes() noexcept {
  return Nodes;
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE constexpr Index
ElementBase<DIM, Derived, Nodes, Faces>::getnFaces() noexcept {
  return Faces;
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE constexpr Index
ElementBase<DIM, Derived, Nodes, Faces>::getnNeighbors() noexcept {
  return Faces;
};

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE Index
ElementBase<DIM, Derived, Nodes, Faces>::getNeighborElems(
    Index iface) const noexcept {
  return neighborElements[iface];
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE void
ElementBase<DIM, Derived, Nodes, Faces>::setNeighborElems(
    Index iface, Index elemIndex) noexcept {
  neighborElements[iface] = elemIndex;
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE Index
ElementBase<DIM, Derived, Nodes, Faces>::getNodes(Index inode) const noexcept {
  return nodes[inode];
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE void ElementBase<DIM, Derived, Nodes, Faces>::setNodes(
    Index inode, Index nodeIndex) noexcept {
  nodes[inode] = nodeIndex;
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE constexpr Index
ElementBase<DIM, Derived, Nodes, Faces>::getMaxNodesFace() const noexcept {
  return crtpCast().getMaxNodesFaceImpl();
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE constexpr Index
ElementBase<DIM, Derived, Nodes, Faces>::getVTKtype() const noexcept {
  return crtpCast().getVTKtypeImpl();
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE real ElementBase<DIM, Derived, Nodes, Faces>::getCoordCenter(
    Index idim) const noexcept {
  return crtpCast().getCoordCenterImpl(idim);
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE void ElementBase<DIM, Derived, Nodes, Faces>::setCoordCenter(
    Index idim, real icoord) noexcept {
  crtpCast().setCoordCenterImpl(idim, icoord);
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE real ElementBase<DIM, Derived, Nodes, Faces>::getCoordFace(
    Index iface, Index idim) const noexcept {
  return crtpCast().getCoordFaceImpl(iface, idim);
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE void ElementBase<DIM, Derived, Nodes, Faces>::setCoordFace(
    Index iface, Index idim, real icoord) noexcept {
  crtpCast().setCoordFaceImpl(iface, idim, icoord);
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE Index
ElementBase<DIM, Derived, Nodes, Faces>::getNodesNeighbor(
    Index inode, Index ineighbor) const noexcept {
  return crtpCast().getNodesNeighborImpl(inode, ineighbor);
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE Index ElementBase<DIM, Derived, Nodes, Faces>::getnNodesFace(
    Index iface) const noexcept {
  return crtpCast().getnNodesFaceImpl(iface);
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE Index
ElementBase<DIM, Derived, Nodes, Faces>::getnNeighborNodes(
    Index inode) const noexcept {
  return crtpCast().getnNeighborNodesImpl(inode);
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE void
ElementBase<DIM, Derived, Nodes, Faces>::changeOrientation() noexcept {
  return crtpCast().changeOrientationImpl();
}

template <int DIM, class Derived, Index Nodes, Index Faces>
RANS_ALWAYS_INLINE Index ElementBase<DIM, Derived, Nodes, Faces>::getFaces(
    Index iface, Index inode) const noexcept {
  return crtpCast().getFacesImpl(iface, inode);
}
}  // namespace preprocess