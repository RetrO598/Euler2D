#pragma once
#include <array>
#include <cassert>
#include <cmath>
#include <common/config.hpp>
namespace preprocess {
template <int DIM>
class Edge {
 public:
  Edge(Index point_0, Index point_1);

  void setCoord(const std::array<std::array<real, DIM>, 2> &icoord);

  real getCoord(Index idim) const noexcept;

  Index getnNodes() const noexcept;

  Index getNode(Index inode) const noexcept;

  std::array<real, DIM> getNormal() const noexcept;

  real getNormal(Index index) const noexcept;

  void getNormal(std::array<real, DIM> &valNormal) const noexcept;

  void setNormal(const std::array<real, DIM> &valNormal);

  void addNormal(const std::array<real, DIM> &valNormal);

  void setZeros();

  real getVolume(const std::array<real, DIM> &coord_edge,
                 const std::array<real, DIM> &coord_face_elem,
                 const std::array<real, DIM> &coord_elem,
                 const std::array<real, DIM> &coord_point) const noexcept;

  real getVolume(const std::array<real, DIM> &coord_edge,
                 const std::array<real, DIM> &coord_elem,
                 const std::array<real, DIM> &coord_point) const noexcept;

  void setNodesCoord(const std::array<real, DIM> &coord_edge,
                     const std::array<real, DIM> &coord_face_elem,
                     const std::array<real, DIM> &coord_elem);

  void setNodesCoord(const std::array<real, DIM> &coord_edge,
                     const std::array<real, DIM> &coord_elem);

 private:
  std::array<real, DIM> coord{};
  std::array<Index, 2> nodes{};
  std::array<real, DIM> normal{};
};

template <int DIM>
Edge<DIM>::Edge(Index point_0, Index point_1) {
  nodes[0] = point_0;
  nodes[1] = point_1;

  for (Index i = 0; i < DIM; ++i) {
    coord[i] = 0.0;
    normal[i] = 0.0;
  }
}

template <int DIM>
void Edge<DIM>::setCoord(const std::array<std::array<real, DIM>, 2> &icoord) {
  for (Index idim = 0; idim < DIM; ++idim) {
    coord[idim] = 0.0;
    for (Index inode = 0; inode < 2; ++inode) {
      coord[idim] += icoord[inode][idim] * 0.5;
    }
  }
}

template <int DIM>
real Edge<DIM>::getCoord(Index idim) const noexcept {
  assert(idim >= 0 && idim < DIM && "idim out of range");
  return coord[idim];
}

template <int DIM>
Index Edge<DIM>::getnNodes() const noexcept {
  return 2;
}

template <int DIM>
Index Edge<DIM>::getNode(Index inode) const noexcept {
  return (inode < 2) ? nodes[inode] : INVALID_INDEX;
}

template <int DIM>
std::array<real, DIM> Edge<DIM>::getNormal() const noexcept {
  return normal;
}

template <int DIM>
real Edge<DIM>::getNormal(Index index) const noexcept {
  return normal[index];
}

template <int DIM>
void Edge<DIM>::getNormal(std::array<real, DIM> &valNormal) const noexcept {
  valNormal = normal;
}

template <int DIM>
void Edge<DIM>::setNormal(const std::array<real, DIM> &valNormal) {
  normal = valNormal;
}

template <int DIM>
void Edge<DIM>::addNormal(const std::array<real, DIM> &valNormal) {
  normal[0] += valNormal[0];
  normal[1] += valNormal[1];
  normal[2] += valNormal[2];
}

template <int DIM>
void Edge<DIM>::setZeros() {
  normal.fill(0.0);
}

template <int DIM>
real Edge<DIM>::getVolume(
    const std::array<real, DIM> &coord_edge,
    const std::array<real, DIM> &coord_face_elem,
    const std::array<real, DIM> &coord_elem,
    const std::array<real, DIM> &coord_point) const noexcept {
  std::array<real, DIM> vec_a, vec_b, vec_c, vec_d;

  for (Index idim = 0; idim < DIM; ++idim) {
    vec_a[idim] = coord_edge[idim] - coord_point[idim];
    vec_b[idim] = coord_face_elem[idim] - coord_point[idim];
    vec_c[idim] = coord_elem[idim] - coord_point[idim];
  }

  vec_d[0] = vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1];
  vec_d[1] = -(vec_a[0] * vec_b[2] - vec_a[2] * vec_b[0]);
  vec_d[2] = vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0];

  real Local_Volume = std::fabs(vec_c[0] * vec_d[0] + vec_c[1] * vec_d[1] +
                                vec_c[2] * vec_d[2]) /
                      6.0;

  return Local_Volume;
}

template <int DIM>
real Edge<DIM>::getVolume(
    const std::array<real, DIM> &coord_edge,
    const std::array<real, DIM> &coord_elem,
    const std::array<real, DIM> &coord_point) const noexcept {
  std::array<real, DIM> vec_a, vec_b;

  for (Index idim = 0; idim < DIM; ++idim) {
    vec_a[idim] = coord_elem[idim] - coord_point[idim];
    vec_b[idim] = coord_edge[idim] - coord_point[idim];
  }

  real Local_Volume =
      0.5 * std::fabs(vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);

  return Local_Volume;
}

template <int DIM>
void Edge<DIM>::setNodesCoord(const std::array<real, DIM> &coord_edge,
                              const std::array<real, DIM> &coord_face_elem,
                              const std::array<real, DIM> &coord_elem) {
  std::array<real, DIM> vec_a, vec_b, dim_normal;

  for (Index idim = 0; idim < DIM; ++idim) {
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
void Edge<DIM>::setNodesCoord(const std::array<real, DIM> &coord_edge,
                              const std::array<real, DIM> &coord_elem) {
  std::array<real, DIM> dim_normal;

  dim_normal[0] = coord_elem[1] - coord_edge[1];
  dim_normal[1] = -(coord_elem[0] - coord_edge[0]);

  normal[0] += dim_normal[0];
  normal[1] += dim_normal[1];
}
}  // namespace preprocess