#pragma once

#include <pre/parameter.h>
#include <yaml-cpp/yaml.h>

#include <string>

namespace preprocess {

class yamlReader {
 public:
  yamlReader(const std::string &filename) : filename(filename) {}

  void read(parameter &param);

 private:
  std::string filename;
};
}  // namespace preprocess