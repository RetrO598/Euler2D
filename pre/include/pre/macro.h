#pragma once

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
  EdgeI *next;
};

struct EdgeList {
  EdgeI *list;
};
} // namespace preprocess