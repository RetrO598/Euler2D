#include <pre/macro.h>
#include <pre/parameter.h>
#include <pre/reader.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
namespace preprocess {

template <typename T>
void extractValues(const std::string& line, const std::regex& re,
                   std::vector<T>& target) {
  for (std::sregex_iterator it(line.begin(), line.end(), re), end; it != end;
       ++it) {
    try {
      if constexpr (std::is_same_v<T, float>) {
        target.push_back(std::stof(it->str()));
      } else if constexpr (std::is_same_v<T, int>) {
        target.push_back(std::stoi(it->str()));
      }
    } catch (const std::exception& e) {
      std::cerr << "Error: Invalid value: " << it->str() << " (" << e.what()
                << ")" << std::endl;
    }
  }
}

std::string simpleReader::readLine() {
  if (file_.eof()) [[unlikely]] {
    return fileEOF;
  }
  last_pos = file_.tellg();
  std::string line;
  std::getline(file_, line);
  return line;
}

std::string simpleReader::readLineFiltered() {
  if (file_.eof()) [[unlikely]] {
    return fileEOF;
  }
  last_pos = file_.tellg();
  std::string line;
  while (std::getline(file_, line)) {
    if (!line.empty() && line.find(commentChar_) == 0) {
      last_pos = file_.tellg();
      continue;
    }
    return line;
  }
  return fileEOF;
}

void simpleReader::customRead(parameter& param) {
  param.title = readLineFiltered();
  param.title.erase(param.title.end() - 1);
  param.gridFile = readLineFiltered();
  param.gridFile.erase(param.gridFile.end() - 1);
  param.flowFieldPlot = readLineFiltered();
  param.surfacePlot = readLineFiltered();
  param.convHistory = readLineFiltered();
  param.restartIn = readLineFiltered();
  param.restartOut = readLineFiltered();

  std::string line = readLineFiltered();

  line = line.substr(0, 3);
  if (line.find("I") != std::string::npos) {
    param.flowtype_ = flowType::Internal;
  } else if (line.find("E") != std::string::npos) {
    param.flowtype_ = flowType::External;
  }

  line = readLineFiltered();
  line = line.substr(0, 3);
  if (line.find("E") != std::string::npos) {
    param.equationtype_ = equationType::Euler;
  } else if (line.find("N") != std::string::npos) {
    param.equationtype_ = equationType::NavierStokes;
  }

  param.gamma = std::stof(readLineFiltered());
  param.Cp = std::stof(readLineFiltered());
  param.Re = std::stof(readLineFiltered());
  param.refVel = std::stof(readLineFiltered());
  param.refRho = std::stof(readLineFiltered());
  param.Prandtl = std::stof(readLineFiltered());

  param.MaInfinity = std::stof(readLineFiltered());
  param.Aoa = std::stof(readLineFiltered());
  param.PsInfinity = std::stof(readLineFiltered());
  param.TsInfinity = std::stof(readLineFiltered());

  param.PtInlet = std::stof(readLineFiltered());
  param.TtInlet = std::stof(readLineFiltered());
  param.flowAngIn = std::stof(readLineFiltered());
  param.PsOutlet = std::stof(readLineFiltered());
  param.approxFlowAngOut = std::stof(readLineFiltered());
  param.PsRatio = std::stof(readLineFiltered());

  param.xRefPoint = std::stof(readLineFiltered());
  param.yRefPoint = std::stof(readLineFiltered());
  param.chord = std::stof(readLineFiltered());

  param.maxIteration = std::stoi(readLineFiltered());
  param.numberIterDump = std::stoi(readLineFiltered());
  param.convTol = std::stof(readLineFiltered());

  line = readLineFiltered();
  line = line.substr(0, 3);
  if (line.find("N") != std::string::npos) {
    param.restart = false;
  } else if (line.find("Y") != std::string::npos) {
    param.restart = true;
  }

  param.CFL = std::stof(readLineFiltered());
  param.imResiSmooth = std::stof(readLineFiltered());
  param.numOfIterSmooth = std::stoi(readLineFiltered());

  line = readLineFiltered();
  line = line.substr(0, 3);
  if (line.find("L") != std::string::npos) {
    param.timestep_ = timeStep::local;
  } else if (line.find("G") != std::string::npos) {
    param.timestep_ = timeStep::global;
  }

  line = readLineFiltered();
  line = line.substr(0, 3);
  if (line.find("N") != std::string::npos) {
    param.preCon = false;
  } else if (line.find("Y") != std::string::npos) {
    param.preCon = true;
  }

  param.preConFactor = std::stof(readLineFiltered());
  param.spatialOrder = std::stoi(readLineFiltered());
  param.limiterCoeff = std::stof(readLineFiltered());
  param.entropyCorreRoe = std::stof(readLineFiltered());

  line = readLineFiltered();
  line = line.substr(0, 3);
  if (line.find("N") != std::string::npos) {
    param.farCorrect = false;
  } else if (line.find("Y") != std::string::npos) {
    param.farCorrect = true;
  }

  param.temperalStages = std::stoi(readLineFiltered());
  param.stageCoeff.reserve(param.temperalStages);
  param.dissipationEval.reserve(param.temperalStages);
  param.dissipationBlend.reserve(param.temperalStages);

  line = readLineFiltered();
  std::regex re(R"(\d+(\.\d+)?)");
  extractValues(line, re, param.stageCoeff);

  line = readLineFiltered();
  extractValues(line, re, param.dissipationBlend);

  line = readLineFiltered();
  extractValues(line, re, param.dissipationEval);

  double gam1 = param.gamma - 1.0;
  double rgas = gam1 * param.Cp / param.gamma;
  if (param.flowtype_ == flowType::External) {
    param.Aoa = param.Aoa * M_PI / 180.0;
    param.RhoInfinity = param.PsInfinity / (rgas * param.TsInfinity);
    param.velInfinity =
        param.MaInfinity * std::sqrt(gam1 * param.Cp * param.TsInfinity);
    param.uInfinity = param.velInfinity * std::cos(param.Aoa);
    param.vInfinity = param.velInfinity * std::sin(param.Aoa);
    param.refMach2 = param.MaInfinity * param.MaInfinity;
    param.refRho = param.RhoInfinity;
    param.refVel = param.velInfinity;
    if (param.equationtype_ == equationType::NavierStokes) {
      param.refVisc =
          param.RhoInfinity * param.velInfinity * param.chord / param.Re;
    } else {
      param.refVisc = 0.0;
    }
  } else {
    double temp = param.TtInlet *
                  std::pow(param.PsOutlet / param.PtInlet, gam1 / param.gamma);
    double mach = std::sqrt(2.0 * ((param.TtInlet / temp) - 1.0) / gam1);

    param.flowAngIn = param.flowAngIn * M_PI / 180.0;
    param.approxFlowAngOut = param.approxFlowAngOut * M_PI / 180.0;
    param.refMach2 = mach * mach;
    if (param.equationtype_ == equationType::NavierStokes) {
      param.refVisc = param.refRho * param.refVel * param.chord / param.Re;
    } else {
      param.refVisc = 0.0;
    }
  }

  param.initialized = true;
}

std::size_t simpleReader::linesCount() const {
  std::ifstream file(filePath, std::ios::in);
  if (!file.is_open()) {
    throw std::ios_base::failure("Failed to open file.");
  }
  size_t lineCount = 0;
  char ch;
  while (file.get(ch)) {
    if (ch == '\n') {
      ++lineCount;
    }
  }
  return lineCount;
}

std::string_view simpleReader::getFilePath() const { return filePath; }
}  // namespace preprocess