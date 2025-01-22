#pragma once

namespace preprocess {
#define fileEOF "end of file"

struct NODE {
  double x, y;
};

struct TRIA {
  int node[3];
};

struct EDGE {
  int nodei, nodej;
};
struct BOUNDFACE {
  int nodei, nodej;
};
struct BOUNDNODE {
  int node, dummy, indexEdge;
};

struct IBOUND {
  int bfaceIndex, bnodeIndex;
};

struct EDGEI {
  int j, edge;
  EDGEI *next;
};

struct EDGELIST {
  EDGEI *list;
};
}  // namespace preprocess