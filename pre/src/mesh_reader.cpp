#include <cstdio>
#include <cstring>
#include <element/hexahedron.hpp>
#include <element/line.hpp>
#include <element/pyramid.hpp>
#include <element/rectangle.hpp>
#include <element/tetrahedron.hpp>
#include <element/triangle.hpp>
#include <element/wedge.hpp>
#include <fstream>
#include <ios>
#include <iostream>
#include <memory>
#include <mesh/mesh_reader.hpp>
#include <sstream>
#include <string>
namespace preprocess {

MeshReader::~MeshReader() {}

inline std::string SU2Reader::trim(const std::string &s) {
  auto start = s.find_first_not_of(" \t\r");
  auto end = s.find_last_not_of(" \t\r");
  return (start == std::string::npos ? "" : s.substr(start, end - start + 1));
}
bool SU2Reader::parseKeyValue(const std::string &line, std::string &key,
                              std::string &value) {
  auto pos = line.find('=');
  if (pos == std::string::npos) return false;

  key = trim(line.substr(0, pos));
  value = trim(line.substr(pos + 1));
  return true;
}

std::unique_ptr<MeshDataBase> SU2Reader::readMesh() const {
  std::ifstream file;
  Index dim = 0;
  file.open(filename, std::ios::in);
  std::string line, key, value;
  std::cout << "Opening file: " << filename << "\n";
  if (!file.is_open()) {
    std::cout << "Mesh file invalid." << "\n";
    return nullptr;
  }

  // First, find the dimension
  while (std::getline(file, line)) {
    if (parseKeyValue(line, key, value)) {
      if (key == "NDIME") {
        dim = std::stoi(value);
        break;  // Found dimension, stop searching
      }
    }
  }

  if (dim == 0) {
    throw std::runtime_error("NDIME not found in mesh file.");
  }

  // Rewind the file to the beginning to parse everything in order
  file.clear();
  file.seekg(0, std::ios::beg);

  if (dim == 2) {
    auto mesh = std::make_unique<MeshData<2>>();
    readMeshDim(file, *mesh);
    return mesh;
  } else if (dim == 3) {
    auto mesh = std::make_unique<MeshData<3>>();
    readMeshDim(file, *mesh);
    return mesh;
  } else {
    throw std::runtime_error("Unsupported dimension: " + std::to_string(dim));
  }

  file.close();
}

template <int DIM>
void SU2Reader::readMeshDim(std::ifstream &file, MeshData<DIM> &mesh) const {
  std::string line, key, value;
  Index nPoints, nElems, index, nMarkers;
  Index vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], vnodes_tetra[4],
      vnodes_hexa[8], vnodes_wedge[6], vnodes_pyramid[5];
  Index nelem_triangle = 0, nelem_quad = 0, nelem_tetra = 0, nelem_hex = 0,
        nelem_wedge = 0, nelem_pyramid = 0;
  Index nelem_edge_bound = 0, nelem_triangle_bound = 0, nelem_quad_bound = 0;
  while (std::getline(file, line)) {
    if (!parseKeyValue(line, key, value)) continue;

    /* --------------------------------------------
     *                NDIME
     * -------------------------------------------- */
    if (key == "NDIME") {
      mesh.nDim = static_cast<Index>(std::stoi(value));
    }
    /* --------------------------------------------
     *                NELEM
     * -------------------------------------------- */
    else if (key == "NELEM") {
      nElems = std::stoi(value);
      mesh.ElemList.reserve(nElems);

      for (Index ielem = 0; ielem < static_cast<Index>(nElems); ++ielem) {
        std::getline(file, line);
        std::istringstream iss(line);

        int elemType;
        iss >> elemType;

        switch (elemType) {
          case TRIANGLE:
            for (Index n = 0; n < TriangleNodes; ++n) {
              iss >> vnodes_triangle[n];
            }
            mesh.ElemList.push_back(
                Triangle<DIM>(static_cast<Index>(vnodes_triangle[0]),
                              static_cast<Index>(vnodes_triangle[1]),
                              static_cast<Index>(vnodes_triangle[2])));
            nelem_triangle++;
            break;
          case RECTANGLE:
            for (Index n = 0; n < RectangleNodes; ++n) {
              iss >> vnodes_quad[n];
            }
            mesh.ElemList.push_back(
                Rectangle<DIM>(static_cast<Index>(vnodes_quad[0]),
                               static_cast<Index>(vnodes_quad[1]),
                               static_cast<Index>(vnodes_quad[2]),
                               static_cast<Index>(vnodes_quad[3])));
            nelem_quad++;
            break;
          case TETRAHEDRON:
            for (Index n = 0; n < TetrahedronNodes; ++n) {
              iss >> vnodes_tetra[n];
            }
            mesh.ElemList.push_back(
                Tetrahedron<DIM>(static_cast<Index>(vnodes_tetra[0]),
                                 static_cast<Index>(vnodes_tetra[1]),
                                 static_cast<Index>(vnodes_tetra[2]),
                                 static_cast<Index>(vnodes_tetra[3])));
            nelem_tetra++;
            break;
          case HEXAHEDRON:
            for (Index n = 0; n < HexahedronNodes; ++n) {
              iss >> vnodes_hexa[n];
            }
            mesh.ElemList.push_back(
                Hexahedron<DIM>(static_cast<Index>(vnodes_hexa[0]),
                                static_cast<Index>(vnodes_hexa[1]),
                                static_cast<Index>(vnodes_hexa[2]),
                                static_cast<Index>(vnodes_hexa[3]),
                                static_cast<Index>(vnodes_hexa[4]),
                                static_cast<Index>(vnodes_hexa[5]),
                                static_cast<Index>(vnodes_hexa[6]),
                                static_cast<Index>(vnodes_hexa[7])));
            nelem_hex++;
            break;
          case WEDGE:
            for (Index n = 0; n < WedgeNodes; ++n) {
              iss >> vnodes_wedge[n];
            }
            mesh.ElemList.push_back(
                Wedge<DIM>(static_cast<Index>(vnodes_wedge[0]),
                           static_cast<Index>(vnodes_wedge[1]),
                           static_cast<Index>(vnodes_wedge[2]),
                           static_cast<Index>(vnodes_wedge[3]),
                           static_cast<Index>(vnodes_wedge[4]),
                           static_cast<Index>(vnodes_wedge[5])));
            nelem_wedge++;
            break;
          case PYRAMID:
            for (Index n = 0; n < PyramidNodes; ++n) {
              iss >> vnodes_pyramid[n];
            }
            mesh.ElemList.push_back(
                Pyramid<DIM>(static_cast<Index>(vnodes_pyramid[0]),
                             static_cast<Index>(vnodes_pyramid[1]),
                             static_cast<Index>(vnodes_pyramid[2]),
                             static_cast<Index>(vnodes_pyramid[3]),
                             static_cast<Index>(vnodes_pyramid[4])));
            nelem_pyramid++;
            break;
          default:
            throw std::runtime_error("Unsupported element type: " +
                                     std::to_string(elemType));
        }
      }
      if (nelem_triangle > 0) {
        std::cout << nelem_triangle << " triangles." << "\n";
      }
      if (nelem_quad > 0) {
        std::cout << nelem_quad << " quadrulaterals." << "\n";
      }
      if (nelem_tetra > 0) {
        std::cout << nelem_tetra << " tetrahedra." << "\n";
      }
      if (nelem_hex > 0) {
        std::cout << nelem_hex << " hexahedra." << "\n";
      }
      if (nelem_wedge > 0) {
        std::cout << nelem_wedge << " prisms." << "\n";
      }
      if (nelem_pyramid > 0) {
        std::cout << nelem_pyramid << " pyramids." << "\n";
      }
    }
    /* --------------------------------------------
     *                NPOIN
     * -------------------------------------------- */
    else if (key == "NPOIN") {
      nPoints = std::stoi(value);
      mesh.PointList.reserve(nPoints);

      for (Index i = 0; i < static_cast<Index>(nPoints); ++i) {
        std::getline(file, line);
        std::istringstream iss(line);

        real x, y, z = 0.0;
        if constexpr (DIM == 2) {
          iss >> x >> y >> index;
          mesh.PointList.emplace_back(x, y, index);
        } else {
          iss >> x >> y >> z >> index;
          mesh.PointList.emplace_back(x, y, z, index);
        }
      }
    }

    /* --------------------------------------------
     *                NMARK
     * -------------------------------------------- */
    else if (key == "NMARK") {
      nMarkers = std::stoi(value);
      mesh.nBound = static_cast<Index>(nMarkers);
      mesh.BoundList.resize(nMarkers);
      mesh.nElemBound.reserve(nMarkers);
      mesh.MarkerTag.reserve(nMarkers);
      mesh.BoundElemIndex.resize(nMarkers);

      for (Index m = 0; m < static_cast<Index>(nMarkers); ++m) {
        // MARKER_TAG
        std::getline(file, line);
        parseKeyValue(line, key, value);
        mesh.MarkerTag.push_back(value);

        // MARKER_ELEMS
        std::getline(file, line);
        parseKeyValue(line, key, value);
        Index nElemB = static_cast<Index>(std::stoi(value));
        mesh.nElemBound.push_back(nElemB);

        mesh.BoundList[m].reserve(nElemB);
        mesh.BoundElemIndex[m].resize(nElemB, INVALID_INDEX);

        // Read marker elements
        for (Index i = 0; i < nElemB; ++i) {
          std::getline(file, line);
          std::istringstream iss(line);

          int elemType;
          iss >> elemType;

          switch (elemType) {
            case LINE:
              for (Index n = 0; n < 2; ++n) {
                iss >> vnodes_edge[n];
              }
              mesh.BoundList[m].push_back(
                  Line<DIM>(static_cast<Index>(vnodes_edge[0]),
                            static_cast<Index>(vnodes_edge[1])));
              nelem_edge_bound++;
              break;
            case TRIANGLE:
              for (Index n = 0; n < TriangleNodes; ++n) {
                iss >> vnodes_triangle[n];
              }
              mesh.BoundList[m].push_back(
                  Triangle<DIM>(static_cast<Index>(vnodes_triangle[0]),
                                static_cast<Index>(vnodes_triangle[1]),
                                static_cast<Index>(vnodes_triangle[2])));
              nelem_triangle_bound++;
              break;
            case RECTANGLE:
              for (Index n = 0; n < RectangleNodes; ++n) {
                iss >> vnodes_quad[n];
              }
              mesh.BoundList[m].push_back(
                  Rectangle<DIM>(static_cast<Index>(vnodes_quad[0]),
                                 static_cast<Index>(vnodes_quad[1]),
                                 static_cast<Index>(vnodes_quad[2]),
                                 static_cast<Index>(vnodes_quad[3])));
              nelem_quad_bound++;
              break;
            default:
              throw std::runtime_error("Unsupported marker element type.");
          }
        }
      }
    }
  }
}

void SU2Reader::writeMesh(const MeshDataBase &mesh,
                          const std::string &outFile) const {
  std::ofstream file;
  file.open(outFile, std::ios::out);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open output file: " + outFile);
  }

  // Set higher buffer size for better I/O performance with large files
  char buffer[65536];
  file.rdbuf()->pubsetbuf(buffer, sizeof(buffer));

  // Lambda: format double in scientific uppercase with 16 decimals using
  // snprintf
  auto fmt = [](real v) -> std::string {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%+.16E", v);
    return std::string(buf);
  };

  // NDIME
  file << "NDIME= " << mesh.nDim << "\n";

  if (mesh.nDim == 2) {
    const auto &mesh2d = dynamic_cast<const MeshData<2> &>(mesh);
    // NELEM (write elements first as requested)
    file << "NELEM= " << static_cast<int>(mesh2d.ElemList.size()) << "\n";
    for (Index ie = 0; ie < mesh2d.ElemList.size(); ++ie) {
      const auto &e = mesh2d.ElemList[ie];
      file << static_cast<int>(e.getVTKtype());
      for (Index inode = 0; inode < e.getnNodes(); ++inode) {
        file << "\t" << static_cast<int>(e.getNodes(inode));
      }
      // append global element index
      file << "\t" << static_cast<int>(ie) << "\n";
    }

    // NPOIN (write points after elements)
    file << "NPOIN= " << static_cast<int>(mesh2d.PointList.size()) << "\n";
    for (Index ip = 0; ip < mesh2d.PointList.size(); ++ip) {
      auto p = mesh2d.PointList[ip];
      file << fmt(p.getCoord(0)) << "\t" << fmt(p.getCoord(1)) << "\t"
           << p.getIndex() << "\n";
    }

    // NMARK (boundary markers)
    file << "NMARK= " << mesh2d.nBound << "\n";
    for (Index m = 0; m < mesh2d.nBound; ++m) {
      file << "MARKER_TAG= " << mesh2d.MarkerTag[m] << "\n";
      file << "MARKER_ELEMS= " << mesh2d.nElemBound[m] << "\n";

      const auto &list = mesh2d.BoundList[m];
      for (const auto &be : list) {
        file << static_cast<int>(be.getVTKtype());
        for (Index inode = 0; inode < be.getnNodes(); ++inode) {
          file << "\t" << static_cast<int>(be.getNodes(inode));
        }
        // For boundary elements do not append global element index
        file << "\n";
      }
    }
  } else if (mesh.nDim == 3) {
    const auto &mesh3d = dynamic_cast<const MeshData<3> &>(mesh);
    // NELEM (write elements first as requested)
    file << "NELEM= " << static_cast<int>(mesh3d.ElemList.size()) << "\n";
    for (Index ie = 0; ie < mesh3d.ElemList.size(); ++ie) {
      const auto &e = mesh3d.ElemList[ie];
      file << static_cast<int>(e.getVTKtype());
      for (Index inode = 0; inode < e.getnNodes(); ++inode) {
        file << "\t" << static_cast<int>(e.getNodes(inode));
      }
      // append global element index
      file << "\t" << static_cast<int>(ie) << "\n";
    }

    // NPOIN (write points after elements)
    file << "NPOIN= " << static_cast<int>(mesh3d.PointList.size()) << "\n";
    for (Index ip = 0; ip < mesh3d.PointList.size(); ++ip) {
      auto p = mesh3d.PointList[ip];
      file << fmt(p.getCoord(0)) << "\t" << fmt(p.getCoord(1)) << "\t"
           << fmt(p.getCoord(2)) << "\t" << p.getIndex() << "\n";
    }

    // NMARK (boundary markers)
    file << "NMARK= " << mesh3d.nBound << "\n";
    for (Index m = 0; m < mesh3d.nBound; ++m) {
      file << "MARKER_TAG= " << mesh3d.MarkerTag[m] << "\n";
      file << "MARKER_ELEMS= " << mesh3d.nElemBound[m] << "\n";

      const auto &list = mesh3d.BoundList[m];
      for (const auto &be : list) {
        file << static_cast<int>(be.getVTKtype());
        for (Index inode = 0; inode < be.getnNodes(); ++inode) {
          file << "\t" << static_cast<int>(be.getNodes(inode));
        }
        // For boundary elements do not append global element index
        file << "\n";
      }
    }
  }

  file.close();
}
}  // namespace preprocess