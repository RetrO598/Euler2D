#include <memory>
#include <pre/geometry.h>
#include <pre/macro.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace preprocess {
Geometry::Geometry(const std::string &filename, const std::string &commentChar)
    : gridReader(filename, commentChar) {
  totNodes = 0;
  phyNodes = 0;
  totEdges = 0;
  phyEdges = 0;
  numTria = 0;
  numBoundSegs = 0;
  numBoundFaces = 0;
  numBoundNodes = 0;
}

void Geometry::ReadGrid() {
  std::string str = gridReader.readLineFiltered();
  std::stringstream(str) >> phyNodes >> numTria >> numBoundSegs;

  BoundTypes.resize(numBoundSegs);
  bname.resize(numBoundSegs);
  ibound.resize(numBoundSegs);

  for (int ib = 0; ib < numBoundSegs; ++ib) {
    str = gridReader.readLineFiltered();
    std::stringstream(str) >> BoundTypes[ib] >> ibound[ib].bfaceIndex >>
        ibound[ib].bnodeIndex;
    ibound[ib].bfaceIndex--;
    ibound[ib].bnodeIndex--;
    bname[ib] = gridReader.readLineFiltered();
  }

  numBoundFaces = ibound[numBoundSegs - 1].bfaceIndex + 1;
  numBoundNodes = ibound[numBoundSegs - 1].bnodeIndex + 1;

  boundaryNode.resize(numBoundNodes);
  boundaryFace.resize(numBoundFaces);

  auto node = BoundaryNode{-777, -777, -777};
  std::fill(boundaryNode.begin(), boundaryNode.end(), node);

  auto face = BoundaryFace{-777, -777};
  std::fill(boundaryFace.begin(), boundaryFace.end(), face);

  int ibegf = 0;
  int iendf = 0;
  int ibegn = 0;
  int iendn = 0;
  for (int i = 0; i < numBoundSegs; ++i) {
    iendf = ibound[i].bfaceIndex;
    iendn = ibound[i].bnodeIndex;
    if (BoundTypes[i] >= 700 && BoundTypes[i] < 800) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        str = gridReader.readLineFiltered();
        std::stringstream(str) >> boundaryNode[ibn].node >>
            boundaryNode[ibn].dummy;
        boundaryNode[ibn].node--;
        boundaryNode[ibn].dummy--;
      }
    } else {
      for (int ibf = ibegf; ibf <= iendf; ++ibf) {
        str = gridReader.readLineFiltered();
        std::stringstream(str) >> boundaryFace[ibf].nodei >>
            boundaryFace[ibf].nodej;
        boundaryFace[ibf].nodei--;
        boundaryFace[ibf].nodej--;
      }
    }
    ibegf = iendf + 1;
    ibegn = iendn + 1;
  }

  for (int i = 0; i < numBoundFaces; ++i) {
    if (boundaryFace[i].nodei < 0 || boundaryFace[i].nodej < 0) {
      throw std::runtime_error("incompletely defined array boundFace[]");
    }
  }

  DummyNodes();

  coords.resize(totNodes);

  for (int i = 0; i < phyNodes; ++i) {
    str = gridReader.readLineFiltered();
    std::stringstream(str) >> coords[i].x >> coords[i].y;
  }

  tria.resize(numTria);

  for (int i = 0; i < numTria; ++i) {
    str = gridReader.readLineFiltered();
    std::stringstream(str) >> tria[i].node[0] >> tria[i].node[1] >>
        tria[i].node[2];
    tria[i].node[0]--;
    tria[i].node[1]--;
    tria[i].node[2]--;
  }

  GenerateEdgeList();

  gridReader.close();
}

void Geometry::DummyNodes() {
  bool flag;
  std::vector<bool> marker;
  marker.resize(phyNodes);

  int ibegf = 0;
  int ibegn = 0;
  int idn = 0;
  int iendf = 0;
  int iendn = 0;
  int itype = 0;
  for (int ib = 0; ib < numBoundSegs; ++ib) {
    iendf = ibound[ib].bfaceIndex;
    iendn = ibound[ib].bnodeIndex;
    itype = BoundTypes[ib];
    flag = (itype >= 100 && itype < 200) || (itype >= 200 && itype < 300) ||
           (itype >= 600 && itype < 700);

    if (itype < 700 || itype >= 800) {
      std::fill(marker.begin(), marker.end(), false);

      for (int i = ibegf; i <= iendf; ++i) {
        marker[boundaryFace[i].nodei] = true;
        marker[boundaryFace[i].nodej] = true;
      }

      for (int i = 0; i < phyNodes; ++i) {
        if (marker[i]) {
          if (ibegn >= numBoundNodes) {
            throw std::runtime_error("Max. number of boundary nodes exceeded.");
          }
          boundaryNode[ibegn].node = i;
          if (flag) {
            boundaryNode[ibegn].dummy = phyNodes + idn;
            boundaryNode[ibegn].indexEdge = -1;
            idn++;
          }
          ibegn++;
        }
      }

      if ((ibegn - 1) != iendn) {
        std::string str = "no. of nodes for boundary " + std::to_string(itype) +
                          " is wrong. It should be " + std::to_string(iendn) +
                          " but it is " + std::to_string(ibegn - 1) + ".";
        throw std::runtime_error(str);
      }
    }

    ibegf = iendf + 1;
    ibegn = iendn + 1;
  }
  marker.clear();
  totNodes = phyNodes + idn;
}

void Geometry::GenerateEdgeList() {
  bool quit;
  int i, j, d, ibn, ie, it, n;
  auto point = std::make_shared<EdgeI>();
  auto prev = std::make_shared<EdgeI>();
  tmpElist.resize(totNodes);

  for (int i = 0; i < phyNodes; ++i) {
    tmpElist[i].list = nullptr;
  }

  phyEdges = 0;

  for (int n = 0; n < 3; ++n) {
    for (int it = 0; it < numTria; ++it) {
      i = tria[it].node[n];
      if (n < 2) {
        j = tria[it].node[n + 1];
      } else {
        j = tria[it].node[0];
      }

      if (i > j) {
        d = i;
        i = j;
        j = d;
      }

      if (tmpElist[i].list == nullptr) {
        phyEdges++;
        tmpElist[i].list = std::make_shared<EdgeI>();
        tmpElist[i].list->j = j;
        tmpElist[i].list->edge = -1;
        tmpElist[i].list->next = nullptr;
      } else if (tmpElist[i].list->j != j) {
        point = tmpElist[i].list->next;
        prev = tmpElist[i].list;
        quit = false;
        do {
          if (point == nullptr) {
            phyEdges++;
            point = std::make_shared<EdgeI>();
            point->j = j;
            point->edge = -1;
            point->next = nullptr;
            prev->next = point;
          }
          if (point->j == j) {
            quit = true;
          }
          prev = point;
          point = point->next;
        } while (!quit);
      }
    }
  }

  totEdges = phyEdges + totNodes - phyNodes;

  edge.resize(totEdges);

  ie = 0;

  for (int i = 0; i < phyNodes; ++i) {
    point = tmpElist[i].list;
    while (point != nullptr) {
      edge[ie].nodei = i;
      edge[ie].nodej = point->j;
      point->edge = ie;
      point = point->next;
      ie++;
    }
  }

  if (ie != phyEdges) {
    throw std::runtime_error("did not get correct number of interior edges.");
  }

  for (ibn = 0; ibn < numBoundNodes; ++ibn) {
    if (boundaryNode[ibn].indexEdge == -1) {
      edge[ie].nodei = boundaryNode[ibn].node;
      edge[ie].nodej = boundaryNode[ibn].dummy;
      boundaryNode[ibn].indexEdge = ie;
      ie++;
    }
  }

  if (ie != totEdges) {
    throw std::runtime_error("did not get the coorect number of dummy edges.");
  }
}

void Geometry::printInfo() {
  std::cout << "No. of interior nodes: " << phyNodes << "\n";
  std::cout << "No. of dummy nodes: " << totNodes - phyNodes << "\n";
  std::cout << "No. of grid cells :" << numTria << "\n";
  std::cout << "No. of interior edges: " << phyEdges << "\n";
  std::cout << "Total number of edges: " << totEdges << "\n";
  std::cout << "No. of boundary faces: " << numBoundFaces << "\n";
  std::cout << "No. of boundary nodes: " << numBoundNodes << "\n";
}

void Geometry::ComputeMetrics() {
  sij.resize(totEdges);
  vol.resize(totNodes);
  sbf.resize(numBoundFaces);
  sproj.resize(totNodes);

  FaceVectorsVolumes();

  tmpElist.clear();

  FaceVectorsVolumesBound();

  CheckMetrics();

  FaceVectorsSymm();

  volumeProjections();
}

void Geometry::FaceVectorsVolumes() {
  int d, i, j, ie, it, n;
  double x1, y1, x2, y2, x3, y3, area, pvol, cx, cy, sx, sy, vprod;
  auto point = std::make_shared<EdgeI>();

  Node node{0.0, 0.0};

  std::fill(sij.begin(), sij.end(), node);
  std::fill(vol.begin(), vol.end(), 0.0);

  for (it = 0; it < numTria; ++it) {
    x1 = coords[tria[it].node[0]].x;
    y1 = coords[tria[it].node[0]].y;
    x2 = coords[tria[it].node[1]].x;
    y2 = coords[tria[it].node[1]].y;
    x3 = coords[tria[it].node[2]].x;
    y3 = coords[tria[it].node[2]].y;

    area = 0.5 * ((x1 - x2) * (y1 + y2) + (x2 - x3) * (y2 + y3) +
                  (x3 - x1) * (y3 + y1));

    pvol = std::fabs(area) / 3.0;

    vol[tria[it].node[0]] += pvol;
    vol[tria[it].node[1]] += pvol;
    vol[tria[it].node[2]] += pvol;

    cx = (x1 + x2 + x3) / 3.0;
    cy = (y1 + y2 + y3) / 3.0;

    for (n = 0; n < 3; ++n) {
      i = tria[it].node[n];
      if (n < 2) {
        j = tria[it].node[n + 1];
      } else {
        j = tria[it].node[0];
      }
      if (i > j) {
        d = i;
        i = j;
        j = d;
      }

      sx = cy - 0.5 * (coords[i].y + coords[j].y);
      sy = -cx + 0.5 * (coords[i].x + coords[j].x);

      vprod =
          sx * (coords[j].x - coords[i].x) + sy * (coords[j].y - coords[i].y);

      if (vprod < 0.0) {
        sx = -sx;
        sy = -sy;
      }

      point = tmpElist[i].list;
      while (point != nullptr) {
        if (point->j == j) {
          ie = point->edge;
          sij[ie].x += sx;
          sij[ie].y += sy;
          break;
        }
        point = point->next;
      }
      if (point == nullptr) {
        throw std::runtime_error("could not find edge to a node.");
      }
    }
  }
}

void Geometry::FaceVectorsVolumesBound() {
  bool flag;
  std::vector<bool> marker;
  int d, i, j, ibegf, ibegn, iendf, iendn, itype, n1, n2, nt1, nt2, nt3;
  int ib, ibf, ibn, ie, it, n;
  std::vector<int> btria;
  double cx, cy, xm, ym, vprod;

  ibegn = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendn = ibound[ib].bnodeIndex;
    if (BoundTypes[ib] >= 700 && BoundTypes[ib] < 800) {
      for (ibn = ibegn; ibn <= iendn; ++ibn) {
        i = boundaryNode[ibn].node;
        j = boundaryNode[ibn].dummy;
        vol[i] += vol[j];
        vol[j] = vol[i];
      }
    }
    ibegn = iendn + 1;
  }

  marker.resize(phyNodes);

  std::fill(marker.begin(), marker.end(), false);

  for (ibf = 0; ibf < numBoundFaces; ++ibf) {
    marker[boundaryFace[ibf].nodei] = true;
    marker[boundaryFace[ibf].nodej] = true;
  }

  btria.resize(numBoundFaces);
  std::fill(btria.begin(), btria.end(), -1);

  for (n = 0; n < 3; ++n) {
    for (it = 0; it < numTria; ++it) {
      i = tria[it].node[n];
      if (n < 2) {
        j = tria[it].node[n + 1];
      } else {
        j = tria[it].node[0];
      }

      if (i > j) {
        d = i;
        i = j;
        j = d;
      }

      if (marker[i] && marker[j]) {
        for (ibf = 0; ibf < numBoundFaces; ++ibf) {
          if (btria[ibf] < 0) {
            n1 = boundaryFace[ibf].nodei;
            n2 = boundaryFace[ibf].nodej;
            if (n1 > n2) {
              d = n1;
              n1 = n2;
              n2 = d;
            }
            if (i == n1 && j == n2) {
              btria[ibf] = it;
            }
          }
        }
      }
    }
  }

  for (ibf = 0; ibf < numBoundFaces; ibf++) {
    if (btria[ibf] < 0)
      throw std::runtime_error(
          "invalid pointer from boundary face to triangle.");
  }

  for (ibf = 0; ibf < numBoundFaces; ++ibf) {
    n1 = boundaryFace[ibf].nodei;
    n2 = boundaryFace[ibf].nodej;
    sbf[ibf].x = coords[n2].y - coords[n1].y;
    sbf[ibf].y = coords[n1].x - coords[n2].x;

    it = btria[ibf];
    nt1 = tria[it].node[0];
    nt2 = tria[it].node[1];
    nt3 = tria[it].node[2];
    cx = (coords[nt1].x + coords[nt2].x + coords[nt3].x) / 3.0;
    cy = (coords[nt1].y + coords[nt2].y + coords[nt3].y) / 3.0;
    xm = 0.5 * (coords[n1].x + coords[n2].x);
    ym = 0.5 * (coords[n1].y + coords[n2].y);
    vprod = sbf[ibf].x * (xm - cx) + sbf[ibf].y * (ym - cy);
    if (vprod < 0.0) {
      sbf[ibf].x = -sbf[ibf].x;
      sbf[ibf].y = -sbf[ibf].y;
    }
  }

  marker.clear();
  btria.clear();

  ibegf = 0;
  ibegn = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendf = ibound[ib].bfaceIndex;
    iendn = ibound[ib].bnodeIndex;
    itype = BoundTypes[ib];
    flag = (itype >= 100 && itype < 200) || (itype >= 200 && itype < 300) ||
           (itype >= 600 && itype < 700);

    if (flag) {
      for (ibf = ibegf; ibf <= iendf; ++ibf) {
        n1 = boundaryFace[ibf].nodei;
        n2 = boundaryFace[ibf].nodej;
        for (ibn = ibegn; ibn <= iendn; ++ibn) {
          if (boundaryNode[ibn].node == n1) {
            ie = boundaryNode[ibn].indexEdge;
            n1 = -777;
            sij[ie].x += 0.5 * sbf[ibf].x;
            sij[ie].y += 0.5 * sbf[ibf].y;
          } else if (boundaryNode[ibn].node == n2) {
            ie = boundaryNode[ibn].indexEdge;
            n2 = -777;
            sij[ie].x += 0.5 * sbf[ibf].x;
            sij[ie].y += 0.5 * sbf[ibf].y;
          }
          if (n1 < 0 && n2 < 0) {
            break;
          }
        }
      }
    }

    ibegf = iendf + 1;
    ibegn = iendn + 1;
  }

  for (ie = phyEdges; ie < totEdges; ++ie) {
    i = edge[ie].nodei;
    j = edge[ie].nodej;
    coords[j].x = coords[i].x;
    coords[j].y = coords[i].y;
    vol[j] = vol[i];
  }
}

void Geometry::CheckMetrics() {
  int i, j, ib, ibf, ibn, ie, ibegf, iendf, ibegn, iendn;
  double volmin, volmax, s, smax;
  std::vector<Node> fvecSum;

  volmin = +1.0e+32;
  volmax = -1.0e+32;

  for (i = 0; i < phyNodes; ++i) {
    volmin = std::min(volmin, vol[i]);
    volmax = std::max(volmax, vol[i]);
  }

  fvecSum.resize(phyNodes);
  Node node{0.0, 0.0};
  std::fill(fvecSum.begin(), fvecSum.end(), node);

  for (ie = 0; ie < phyEdges; ++ie) {
    i = edge[ie].nodei;
    j = edge[ie].nodej;
    fvecSum[i].x += sij[ie].x;
    fvecSum[i].y += sij[ie].y;
    fvecSum[j].x -= sij[ie].x;
    fvecSum[j].y -= sij[ie].y;
  }

  ibegf = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendf = ibound[ib].bfaceIndex;
    if (BoundTypes[ib] < 700 || BoundTypes[ib] >= 800) {
      for (ibf = ibegf; ibf <= iendf; ++ibf) {
        i = boundaryFace[ibf].nodei;
        j = boundaryFace[ibf].nodej;
        fvecSum[i].x += 0.5 * sbf[ibf].x;
        fvecSum[i].y += 0.5 * sbf[ibf].y;
        fvecSum[j].x += 0.5 * sbf[ibf].x;
        fvecSum[j].y += 0.5 * sbf[ibf].y;
      }
    }

    ibegf = iendf + 1;
  }

  ibegn = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendn = ibound[ib].bnodeIndex;
    if (BoundTypes[ib] >= 700 && BoundTypes[ib] < 800) {
      for (ibn = ibegn; ibn <= iendn; ++ibn) {
        i = boundaryNode[ibn].node;
        j = boundaryNode[ibn].dummy;
        fvecSum[i].x += fvecSum[j].x;
        fvecSum[i].y += fvecSum[j].y;
        fvecSum[j].x = fvecSum[i].x;
        fvecSum[j].y = fvecSum[i].y;
      }
    }
    ibegn = iendn + 1;
  }

  smax = -1.0e+32;
  for (i = 0; i < phyNodes; ++i) {
    s = std::sqrt(fvecSum[i].x * fvecSum[i].x + fvecSum[i].y * fvecSum[i].y);
    smax = std::max(smax, s);
  }

  std::cout << " max. sum(S) = " << smax << "\n";
  std::cout << " min. volume = " << volmin << "\n";
  std::cout << " max. volume = " << volmax << "\n";
}

void Geometry::FaceVectorsSymm() {
  int i, j, ib, ibf, ibn, ie, ibegf, ibegn, iendf, iendn;
  std::vector<int> marker;
  double sx, sy;

  marker.resize(phyNodes);

  std::fill(marker.begin(), marker.end(), -1);

  ibegf = 0;
  ibegn = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendf = ibound[ib].bfaceIndex;
    iendn = ibound[ib].bnodeIndex;
    if (BoundTypes[ib] >= 500 && BoundTypes[ib] < 600) {
      sx = 0.0;
      sy = 0.0;
      for (ibf = ibegf; ibf <= iendf; ++ibf) {
        sx += sbf[ibf].x;
        sy += sbf[ibf].y;
      }
      if (std::abs(sx) > std::abs(sy)) {
        BoundTypes[ib] = 501;
      } else {
        BoundTypes[ib] = 502;
      }

      for (ibn = ibegn; ibn <= iendn; ++ibn) {
        marker[boundaryNode[ibn].node] = BoundTypes[ib] - 500;
      }
    }
    ibegf = iendf + 1;
    ibegn = iendn + 1;
  }

  for (ie = 0; ie < phyEdges; ++ie) {
    i = edge[ie].nodei;
    j = edge[ie].nodej;
    if (marker[i] != -1 && marker[j] != -1) {
      if (marker[i] < 2) // x=const. plane
        sij[ie].x = 0.0;
      else // y=const. plane
        sij[ie].y = 0.0;
    }
  }
}

void Geometry::volumeProjections() {
  int ibegn, iendn, ibegf, iendf;
  int i, j, ib, ibf, ibn, ie;
  double sx, sy;

  for (i = 0; i < totNodes; ++i) {
    sproj[i].x = 0.0;
    sproj[i].y = 0.0;
  }

  for (ie = 0; ie < phyEdges; ++ie) {
    i = edge[ie].nodei;
    j = edge[ie].nodej;
    sx = 0.5 * std::fabs(sij[ie].x);
    sy = 0.5 * std::fabs(sij[ie].y);
    sproj[i].x += sx;
    sproj[i].y += sy;
    sproj[j].x += sx;
    sproj[j].y += sy;
  }

  ibegf = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendf = ibound[ib].bfaceIndex;
    if (BoundTypes[ib] < 700 || BoundTypes[ib] >= 800) {
      for (ibf = ibegf; ibf <= iendf; ++ibf) {
        i = boundaryFace[ibf].nodei;
        j = boundaryFace[ibf].nodej;
        sx = 0.25 * std::fabs(sbf[ibf].x);
        sy = 0.25 * std::fabs(sbf[ibf].y);
        sproj[i].x += sx;
        sproj[i].y += sy;
        sproj[j].x += sx;
        sproj[j].y += sy;
      }
    }
    ibegf = iendf + 1;
  }

  ibegn = 0;
  for (ib = 0; ib < numBoundSegs; ++ib) {
    iendn = ibound[ib].bnodeIndex;
    if (BoundTypes[ib] >= 700 && BoundTypes[ib] < 800) {
      for (ibn = ibegn; ibn <= iendn; ibn++) {
        i = boundaryNode[ibn].node;
        j = boundaryNode[ibn].dummy;
        sproj[i].x += sproj[j].x;
        sproj[i].y += sproj[j].y;
        sproj[j].x = sproj[i].x;
        sproj[j].y = sproj[i].y;
      }
    }
    ibegn = iendn + 1;
  }
}

void Geometry::outputMeshInfo() {
  std::ofstream outputFile1("sij.txt");

  for (int i = 0; i < totEdges; ++i) {
    outputFile1 << sij[i].x << " " << sij[i].y << "\n";
  }

  outputFile1.close();

  std::ofstream outputFile2("edge.txt");
  for (int i = 0; i < totEdges; ++i) {
    outputFile2 << edge[i].nodei << " " << edge[i].nodej << "\n";
  }

  outputFile2.close();

  std::ofstream outputFile3("vol.txt");
  for (int i = 0; i < totNodes; ++i) {
    outputFile3 << vol[i] << "\n";
  }

  outputFile3.close();

  std::ofstream outputFile4("sbf.txt");
  for (int i = 0; i < numBoundFaces; ++i) {
    outputFile4 << sbf[i].x << " " << sbf[i].y << "\n";
  }

  outputFile4.close();

  std::ofstream outputFile5("sproj.txt");
  for (int i = 0; i < totNodes; ++i) {
    outputFile5 << sproj[i].x << " " << sproj[i].y << "\n";
  }

  outputFile5.close();
}
} // namespace preprocess