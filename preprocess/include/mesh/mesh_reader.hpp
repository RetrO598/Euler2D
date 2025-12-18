#pragma once

#include <sys/stat.h>

#include <memory>
#include <mesh/mesh_data.hpp>
#include <string>
namespace preprocess {
class MeshReader {
 public:
  virtual std::unique_ptr<MeshDataBase> readMesh() const = 0;
  virtual ~MeshReader() = 0;
};

class SU2Reader : public MeshReader {
 public:
  SU2Reader(std::string meshfile) : filename(meshfile) {}
  std::unique_ptr<MeshDataBase> readMesh() const override;

  void writeMesh(const MeshDataBase &mesh, const std::string &outFile) const;

  static std::string trim(const std::string &s);

  static bool parseKeyValue(const std::string &line, std::string &key,
                            std::string &value);

  ~SU2Reader() = default;

 private:
  std::string filename;
  template <int DIM>
  void readMeshDim(std::ifstream &file, MeshData<DIM> &mesh) const;
};
}  // namespace preprocess