#pragma once
#include <memory>
namespace preprocess {
#define fileEOF "end of file"

struct Node {
  double x, y;
};

struct Tria {
  int node[3];
};

struct Edge {
  int nodei, nodej;
};
struct BoundaryFace {
  int nodei, nodej;
};
struct BoundaryNode {
  int node, dummy, indexEdge;
};

struct idBoundary {
  int bfaceIndex, bnodeIndex;
};

struct EdgeI {
  int j, edge;
  std::shared_ptr<EdgeI> next;
};

struct EdgeList {
  std::shared_ptr<EdgeI> list;
};
} // namespace preprocess