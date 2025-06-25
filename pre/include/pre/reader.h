#pragma once

#include <pre/parameter.h>
#include <yaml-cpp/yaml.h>

#include <cstddef>
#include <fstream>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <string_view>

namespace preprocess {

class simpleReader {
 public:
  simpleReader(const std::string &filename,
               const std::string &commentChar = "#")
      : filePath(filename),
        file_(filename),
        last_pos(0),
        commentChar_(commentChar) {
    if (!file_.is_open()) {
      throw std::runtime_error("Unable to open file: " + filePath);
    }
  };

  simpleReader() = delete;

  ~simpleReader() {
    if (file_.is_open()) {
      file_.close();
    }
  }

  inline void close() {
    if (file_.is_open()) {
      file_.close();
    }
  }

  std::string readLine();

  std::string readLineFiltered();

  virtual void customRead(parameter &param);

  std::size_t linesCount() const;

  std::string_view getFilePath() const;

 private:
  std::string filePath;
  std::ifstream file_;
  std::streampos last_pos;
  std::string commentChar_;
};

class yamlReader {
 public:
  yamlReader(const std::string &filename) : filename(filename) {}

  void read(parameter &param);

 private:
  std::string filename;
};
}  // namespace preprocess