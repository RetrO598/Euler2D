#pragma once

#include <array>
#include <common/config.hpp>

namespace preprocess {
template <int DIM>
class Vertex {
 public:
  Vertex(Index ipoint);

  void setNodesCoord(const std::array<real, DIM> &coord_edge,
                     const std::array<real, DIM> &coord_face_elem,
                     const std::array<real, DIM> &coord_elem);

  void setNodesCoord(const std::array<real, DIM> &coord_edge,
                     const std::array<real, DIM> &coord_elem);

  void addNormal(const std::array<real, DIM> &face_normal);

  Index getnNodes() const noexcept;

  Index getNode() const noexcept;

  std::array<real, DIM> getNormal() const noexcept;

  void getNormal(std::array<real, DIM> &valNormal) const noexcept;

  real getNormal(Index index) const noexcept;

  void setNormal(const std::array<real, DIM> &valNormal);

  std::array<real, DIM> getCoord() const noexcept;

  real getCoord(int idim) const noexcept;

  void setCoord(const std::array<real, DIM> &valCoord);

  Index getNormalNeighbor() const noexcept;

  void setNormalNeighbor(Index pointNormal) noexcept;

  Index getPeriodicPair() const noexcept;

  void setPeriodicPair(Index nodeIndex) noexcept;

 private:
  Index node;
  std::array<real, DIM> normal{};
  Index normalNeighbor;
  Index PeriodicPair;
  std::array<real, DIM> coord{};
};

template <int DIM>
Vertex<DIM>::Vertex(Index ipoint)
    : node(ipoint),
      normalNeighbor(INVALID_INDEX),
      PeriodicPair(INVALID_INDEX) {}

template <int DIM>
void Vertex<DIM>::setNodesCoord(const std::array<real, DIM> &coord_edge,
                                const std::array<real, DIM> &coord_face_elem,
                                const std::array<real, DIM> &coord_elem) {
  std::array<real, DIM> vec_a, vec_b, dim_normal;

  for (int idim = 0; idim < DIM; ++idim) {
    vec_a[idim] = coord_elem[idim] - coord_edge[idim];
    vec_b[idim] = coord_face_elem[idim] - coord_edge[idim];
  }

  dim_normal[0] = 0.5 * (vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1]);
  dim_normal[1] = -0.5 * (vec_a[0] * vec_b[2] - vec_a[2] * vec_b[0]);
  dim_normal[2] = 0.5 * (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);

  normal[0] += dim_normal[0];
  normal[1] += dim_normal[1];
  normal[2] += dim_normal[2];
}

template <int DIM>
void Vertex<DIM>::setNodesCoord(const std::array<real, DIM> &coord_edge,
                                const std::array<real, DIM> &coord_elem) {
  std::array<real, DIM> dim_normal;

  dim_normal[0] = coord_elem[1] - coord_edge[1];
  dim_normal[1] = -(coord_elem[0] - coord_edge[0]);

  normal[0] += dim_normal[0];
  normal[1] += dim_normal[1];
}

template <int DIM>
void Vertex<DIM>::addNormal(const std::array<real, DIM> &face_normal) {
  normal[0] += face_normal[0];
  normal[1] += face_normal[1];
  normal[2] += face_normal[2];
}

template <int DIM>
Index Vertex<DIM>::getnNodes() const noexcept {
  return 1;
}

template <int DIM>
Index Vertex<DIM>::getNode() const noexcept {
  return node;
}

template <int DIM>
std::array<real, DIM> Vertex<DIM>::getNormal() const noexcept {
  return normal;
}

template <int DIM>
real Vertex<DIM>::getNormal(Index index) const noexcept {
  return normal[index];
}

template <int DIM>
void Vertex<DIM>::getNormal(std::array<real, DIM> &valNormal) const noexcept {
  valNormal = normal;
}

template <int DIM>
void Vertex<DIM>::setNormal(const std::array<real, DIM> &valNormal) {
  normal = valNormal;
}

template <int DIM>
std::array<real, DIM> Vertex<DIM>::getCoord() const noexcept {
  return coord;
}

template <int DIM>
real Vertex<DIM>::getCoord(int idim) const noexcept {
  if (idim < 0 || idim >= DIM) {
    return 0.0;
  }
  return coord[idim];
}

template <int DIM>
void Vertex<DIM>::setCoord(const std::array<real, DIM> &valCoord) {
  coord = valCoord;
}

template <int DIM>
Index Vertex<DIM>::getNormalNeighbor() const noexcept {
  return normalNeighbor;
}

template <int DIM>
void Vertex<DIM>::setNormalNeighbor(Index pointNormal) noexcept {
  normalNeighbor = pointNormal;
}

template <int DIM>
Index Vertex<DIM>::getPeriodicPair() const noexcept {
  return PeriodicPair;
}

template <int DIM>
void Vertex<DIM>::setPeriodicPair(Index nodeIndex) noexcept {
  PeriodicPair = nodeIndex;
}
}  // namespace preprocess