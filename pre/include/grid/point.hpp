#pragma once

#include <array>
#include <cassert>
#include <common/config.hpp>
#include <vector>

namespace preprocess {
template <int DIM>
class Point {
 public:
  template <int D = DIM, typename = std::enable_if_t<D == 2>>
  Point(real coord_0, real coord_1, Index index);
  template <int D = DIM, typename = std::enable_if_t<D == 3>>
  Point(real coord_0, real coord_1, real coord_2, Index index);

  Index getIndex() const noexcept;

  void setIndex(Index val) noexcept;

  void setPoint(Index val);

  Index getPoint(Index ipoint) const noexcept;

  void initBoundary(Index val);

  void setElem(Index val);

  void resetBoundary();

  void resetElem();

  void resetPoint();

  real getCoord(Index idim) const noexcept;

  std::array<real, DIM> getCoord() const noexcept;

  void setCoord(Index idim, real icoord);

  void addCoord(Index idim, real icoord);

  void setCoord(const std::array<real, DIM> &coord);

  void setnElem(Index nelem);

  Index getnElem() const noexcept;

  void setEdge(Index iedge, Index valEdge);

  Index getEdge(Index iedge) const noexcept;

  Index getElem(Index ielem) const noexcept;

  void setnPoint(Index npoint);

  Index getnPoint() const noexcept;

  real getVolume() const noexcept;

  bool getBoundary() const noexcept;

  void setBoundary(Index val);

  void setPhysicalBoundary(bool val);

  bool getPhysicalBoundary() const noexcept;

  void setSolidBoundary(bool val);

  bool getSolidBoundary() const noexcept;

  void setDomain(bool val) noexcept;

  bool getDomain() const noexcept;

  void addVolume(real val);

  void setVolume(real val);

  void setVertex(Index iVertex, Index valVertex);

  Index getVertex(Index iVertex) const noexcept;

 private:
  Index index;
  Index nElem;
  Index nPoint;
  std::vector<Index> elems;
  std::vector<Index> points;
  std::vector<Index> edges;
  real volume;
  bool domain, boundary, physicalBoundary, solidBoundary;
  std::vector<Index> vertex;
  std::array<real, DIM> coord{};
};

template <int DIM>
template <int D, typename>
Point<DIM>::Point(real coord_0, real coord_1, Index index)
    : index(int_to_index(index)),
      nElem(0),
      nPoint(0),
      volume(0.0),
      domain(true),
      boundary(false),
      physicalBoundary(false),
      solidBoundary(false),
      coord{coord_0, coord_1} {}

template <int DIM>
template <int D, typename>
Point<DIM>::Point(real coord_0, real coord_1, real coord_2, Index index)
    : index(int_to_index(index)),
      nElem(0),
      nPoint(0),
      volume(0.0),
      domain(true),
      boundary(false),
      physicalBoundary(false),
      solidBoundary(false),
      coord{coord_0, coord_1, coord_2} {}

template <int DIM>
Index Point<DIM>::getIndex() const noexcept {
  return index;
}

template <int DIM>
void Point<DIM>::setIndex(Index val) noexcept {
  index = val;
}

template <int DIM>
void Point<DIM>::setPoint(Index val) {
  bool new_point = true;
  for (Index iPoint = 0; iPoint < points.size(); ++iPoint) {
    if (points[iPoint] == val) {
      new_point = false;
      break;
    }
  }

  if (new_point) {
    points.push_back(val);
    // keep a parallel slot in `edges` for edge/index bookkeeping
    edges.push_back(INVALID_INDEX);
    nPoint = points.size();
  }
}

template <int DIM>
Index Point<DIM>::getPoint(Index ipoint) const noexcept {
  if (ipoint < points.size()) {
    return points[ipoint];
  }
  return INVALID_INDEX;
}

template <int DIM>
void Point<DIM>::initBoundary(Index val) {
  if (!boundary) {
    vertex.resize(val);
    for (Index i = 0; i < val; ++i) {
      vertex[i] = INVALID_INDEX;
    }
  }

  boundary = true;
}

template <int DIM>
void Point<DIM>::setElem(Index val) {
  elems.push_back(val);
  nElem = elems.size();
}

template <int DIM>
void Point<DIM>::resetBoundary() {
  if (!vertex.empty()) {
    vertex.clear();
  }
  boundary = false;
}

template <int DIM>
void Point<DIM>::resetElem() {
  elems.clear();
  nElem = 0;
}

template <int DIM>
void Point<DIM>::resetPoint() {
  points.clear();
  edges.clear();
  nPoint = 0;
}

template <int DIM>
real Point<DIM>::getCoord(Index idim) const noexcept {
  assert(idim >= 0 && idim < DIM && "idim out of range");
  return coord[idim];
}

template <int DIM>
std::array<real, DIM> Point<DIM>::getCoord() const noexcept {
  return coord;
}

template <int DIM>
void Point<DIM>::setCoord(Index idim, real icoord) {
  assert(idim >= 0 && idim < DIM && "idim out of range");
  coord[idim] = icoord;
}

template <int DIM>
void Point<DIM>::addCoord(Index idim, real icoord) {
  assert(idim >= 0 && idim < DIM && "idim out of range");
  coord[idim] += icoord;
}

template <int DIM>
void Point<DIM>::setCoord(const std::array<real, DIM> &icoord) {
  coord = icoord;
}

template <int DIM>
void Point<DIM>::setnElem(Index nelem) {
  nElem = nelem;
}

template <int DIM>
Index Point<DIM>::getnElem() const noexcept {
  return nElem;
}

template <int DIM>
void Point<DIM>::setEdge(Index iedge, Index valEdge) {
  if (iedge < edges.size()) {
    edges[iedge] = valEdge;
  }
}

template <int DIM>
Index Point<DIM>::getEdge(Index iedge) const noexcept {
  if (iedge < edges.size()) {
    return edges[iedge];
  }
  return INVALID_INDEX;
}

template <int DIM>
Index Point<DIM>::getElem(Index ielem) const noexcept {
  if (ielem < elems.size()) {
    return elems[ielem];
  }
  return INVALID_INDEX;
}

template <int DIM>
void Point<DIM>::setnPoint(Index npoint) {
  nPoint = npoint;
}

template <int DIM>
Index Point<DIM>::getnPoint() const noexcept {
  return nPoint;
}

template <int DIM>
real Point<DIM>::getVolume() const noexcept {
  return volume;
}

template <int DIM>
bool Point<DIM>::getBoundary() const noexcept {
  return boundary;
}

template <int DIM>
void Point<DIM>::setBoundary(Index val) {
  if (!boundary) {
    vertex.resize(val);
    for (Index i = 0; i < val; ++i) {
      vertex[i] = INVALID_INDEX;
    }
  }
  boundary = true;
}

template <int DIM>
void Point<DIM>::setPhysicalBoundary(bool val) {
  physicalBoundary = val;
}

template <int DIM>
bool Point<DIM>::getPhysicalBoundary() const noexcept {
  return physicalBoundary;
}

template <int DIM>
void Point<DIM>::setSolidBoundary(bool val) {
  solidBoundary = val;
}

template <int DIM>
bool Point<DIM>::getSolidBoundary() const noexcept {
  return solidBoundary;
}

template <int DIM>
void Point<DIM>::setDomain(bool val) noexcept {
  domain = val;
}

template <int DIM>
bool Point<DIM>::getDomain() const noexcept {
  return domain;
}

template <int DIM>
void Point<DIM>::addVolume(real val) {
  volume += val;
}

template <int DIM>
void Point<DIM>::setVolume(real val) {
  volume = val;
}

template <int DIM>
void Point<DIM>::setVertex(Index iVertex, Index valVertex) {
  if (iVertex < vertex.size() && boundary) {
    vertex[iVertex] = valVertex;
  }
}

template <int DIM>
Index Point<DIM>::getVertex(Index iVertex) const noexcept {
  if (iVertex < vertex.size() && boundary) {
    return vertex[iVertex];
  }
  return INVALID_INDEX;
}
}  // namespace preprocess