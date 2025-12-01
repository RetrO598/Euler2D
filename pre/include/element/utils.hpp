#pragma once

#include <cstddef>
namespace preprocess {

constexpr static std::size_t RectangleNodes = 4;
constexpr static std::size_t RectangleFaces = 4;

constexpr static std::size_t TriangleNodes = 3;
constexpr static std::size_t TriangleFaces = 3;

constexpr static std::size_t LineNodes = 2;
constexpr static std::size_t LineFaces = 1;

constexpr static std::size_t TetrahedronNodes = 4;
constexpr static std::size_t TetrahedronFaces = 4;

constexpr static std::size_t HexahedronNodes = 8;
constexpr static std::size_t HexahedronFaces = 6;

constexpr static std::size_t WedgeNodes = 6;
constexpr static std::size_t WedgeFaces = 5;

constexpr static std::size_t PyramidNodes = 5;
constexpr static std::size_t PyramidFaces = 5;

enum ELEMENTTYPE {
  VERTEX = 1,
  LINE = 3,
  TRIANGLE = 5,
  RECTANGLE = 9,
  TETRAHEDRON = 10,
  HEXAHEDRON = 12,
  WEDGE = 13,
  PYRAMID = 14
};
}  // namespace preprocess