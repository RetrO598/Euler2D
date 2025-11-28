#include <pre/macro.h>
#include <pre/parameter.h>
#include <pre/reader.h>
#include <yaml-cpp/node/parse.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
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

void yamlReader::read(parameter& param) {
  YAML::Node config = YAML::LoadFile(filename);
  auto get = [&](const std::string& key) -> YAML::Node { return config[key]; };

  // 文件与输出
  param.title = get("title").as<std::string>();
  param.gridFile = get("files")["grid_file"].as<std::string>();
  param.flowFieldPlot = get("files")["flow_field"].as<std::string>();
  param.surfacePlot = get("files")["surface_quantities"].as<std::string>();
  param.convHistory = get("files")["convergence_history"].as<std::string>();
  param.restartIn = get("files")["restart_in"].as<std::string>();
  param.restartOut = get("files")["restart_out"].as<std::string>();

  // 物理设置
  std::string flow = get("physics")["general"]["flow_type"].as<std::string>();
  param.flowtype_ = (flow == "I" ? flowType::Internal : flowType::External);

  std::string eqn =
      get("physics")["general"]["equation_type"].as<std::string>();
  param.equationtype_ =
      (eqn == "E" ? equationType::Euler : equationType::NavierStokes);

  param.gamma = get("physics")["general"]["gamma"].as<double>();
  param.Cp = get("physics")["general"]["cp"].as<double>();
  param.Re = get("physics")["general"]["reynolds_number"].as<double>();
  param.refVel = get("physics")["general"]["reference_velocity"].as<double>();
  param.refRho = get("physics")["general"]["reference_density"].as<double>();
  param.Prandtl = get("physics")["general"]["prandtl_number"].as<double>();

  // 外部流参数
  param.MaInfinity = get("physics")["external_flow"]["mach"].as<double>();
  param.Aoa = get("physics")["external_flow"]["angle_of_attack"].as<double>();
  param.PsInfinity =
      get("physics")["external_flow"]["static_pressure"].as<double>();
  param.TsInfinity =
      get("physics")["external_flow"]["static_temperature"].as<double>();

  // 内部流参数
  param.PtInlet =
      get("physics")["internal_flow"]["total_pressure_inlet"].as<double>();
  param.TtInlet =
      get("physics")["internal_flow"]["total_temperature_inlet"].as<double>();
  param.flowAngIn =
      get("physics")["internal_flow"]["flow_angle_inlet"].as<double>();
  param.PsOutlet =
      get("physics")["internal_flow"]["static_pressure_outlet"].as<double>();
  param.approxFlowAngOut =
      get("physics")["internal_flow"]["flow_angle_outlet"].as<double>();
  param.PsRatio = get("physics")["internal_flow"]["pressure_ratio_inlet_outlet"]
                      .as<double>();

  auto inflow =
      get("boundary_condition")["Inflow"].as<std::vector<std::string>>();
  if (!inflow.empty()) {
    for (auto i : inflow) {
      param.boundaryMap.insert({i, BoundaryType::Inflow});
    }
  }

  auto outflow =
      get("boundary_condition")["Outflow"].as<std::vector<std::string>>();
  if (!outflow.empty()) {
    for (auto i : outflow) {
      param.boundaryMap.insert({i, BoundaryType::Outflow});
    }
  }

  auto farfield =
      get("boundary_condition")["Farfield"].as<std::vector<std::string>>();
  if (!farfield.empty()) {
    for (auto i : farfield) {
      param.boundaryMap.insert({i, BoundaryType::Farfield});
    }
  }

  auto eulerWall =
      get("boundary_condition")["EulerWall"].as<std::vector<std::string>>();
  if (!eulerWall.empty()) {
    for (auto i : eulerWall) {
      param.boundaryMap.insert({i, BoundaryType::EulerWall});
    }
  }

  auto noslipWall =
      get("boundary_condition")["NoSlipWall"].as<std::vector<std::string>>();
  if (!noslipWall.empty()) {
    for (auto i : noslipWall) {
      param.boundaryMap.insert({i, BoundaryType::NoSlipWall});
    }
  }

  auto periodic =
      get("boundary_condition")["Periodic"].as<std::vector<std::string>>();
  if (!periodic.empty()) {
    for (auto i : periodic) {
      param.boundaryMap.insert({i, BoundaryType::Periodic});
    }
  }

  auto symmetry =
      get("boundary_condition")["Symmetric"].as<std::vector<std::string>>();
  if (!symmetry.empty()) {
    for (auto i : symmetry) {
      param.boundaryMap.insert({i, BoundaryType::Symmetric});
    }
  }

  for (auto i : param.boundaryMap) {
    std::cout << "Boundary: " << i.first << " with type: " << int(i.second)
              << "\n";
  }

  // 几何参数
  param.xRefPoint = get("geometry_reference")["x_ref"].as<double>();
  param.yRefPoint = get("geometry_reference")["y_ref"].as<double>();
  param.chord = get("geometry_reference")["chord_length"].as<double>();

  // 迭代控制
  param.maxIteration = get("iteration_control")["max_iterations"].as<int>();
  param.numberIterDump = get("iteration_control")["dump_frequency"].as<int>();
  param.convTol = get("iteration_control")["tolerance"].as<double>();
  param.restart = get("iteration_control")["restart"].as<bool>();

  // 数值控制
  auto num = get("numerical_parameters");
  param.CFL = num["cfl"].as<double>();
  param.imResiSmooth = num["smoothing_coeff"].as<double>();
  param.numOfIterSmooth = num["smoothing_jacobi_iters"].as<int>();

  std::string tstep = num["time_stepping"].as<std::string>();
  param.timestep_ = (tstep == "L" ? timeStep::local : timeStep::global);

  param.preCon = num["low_mach_preconditioning"].as<bool>();
  param.preConFactor = num["preconditioning_K"].as<double>();
  param.spatialOrder = num["roe_order"].as<int>();
  param.limiterCoeff = num["limiter_coefficient"].as<double>();
  param.entropyCorreRoe = num["entropy_fix_coefficient"].as<double>();
  param.farCorrect = num["vortex_farfield_correction"].as<bool>();

  // 多阶段时间推进
  param.temperalStages = num["stages"].as<int>();
  param.stageCoeff = num["stage_coefficients"].as<std::vector<double>>();
  param.dissipationBlend =
      num["dissipation_blending"].as<std::vector<double>>();
  param.dissipationEval = num["dissipation_evaluation"].as<std::vector<int>>();

  // 派生量计算
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
    param.refVisc =
        (param.equationtype_ == equationType::NavierStokes)
            ? param.RhoInfinity * param.velInfinity * param.chord / param.Re
            : 0.0;
  } else {
    double temp = param.TtInlet *
                  std::pow(param.PsOutlet / param.PtInlet, gam1 / param.gamma);
    double mach = std::sqrt(2.0 * ((param.TtInlet / temp) - 1.0) / gam1);
    param.flowAngIn *= M_PI / 180.0;
    param.approxFlowAngOut *= M_PI / 180.0;
    param.refMach2 = mach * mach;
    param.refVisc = (param.equationtype_ == equationType::NavierStokes)
                        ? param.refRho * param.refVel * param.chord / param.Re
                        : 0.0;
  }

  param.initialized = true;
}
}  // namespace preprocess