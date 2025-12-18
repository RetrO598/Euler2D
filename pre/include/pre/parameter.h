#pragma once
#include <array>
#include <cstddef>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace preprocess {

enum class BoundaryType {
  Inflow,
  Outflow,
  Farfield,
  EulerWall,
  NoSlipWall,
  Periodic,
  Symmetric
};

enum class ConvectionScheme { ROE, SLAU2, AUSM, AUSMUP2 };

enum class TemporalScheme { RungeKutta, LUSGS };

enum class Limiter { VenkataKrishnan, NishikawaR3 };

enum class flowType { Internal, External };

enum class equationType { Euler, NavierStokes, RANS };

enum class timeStep { local, global };

struct PeriodicInfo {
  std::string master_name;
  std::string slave_name;
  std::array<double, 2> center;
  double rotate;
  std::array<double, 2> translate;
};

struct parameter {
  // explicit parameter(){};

  // parameter(const parameter &other) = delete;

  // parameter &operator=(const parameter &other) = delete;

  // parameter(parameter &&other) = delete;

  // parameter &operator=(parameter &&other) = delete;

  // ~parameter() = default;

  void printParameters();

  bool initialized = false;

  // file name
  std::string title, gridFile, flowFieldPlot, surfacePlot, convHistory,
      restartIn, restartOut;

  // E = external flow, I = internal flow
  flowType flowtype_;

  // E = Euler, N = Navier-Stokes
  equationType equationtype_;

  // some parameters
  double gamma, Cp, Re, refVel, refRho, Prandtl, refVisc, refMach2;

  // Physics parameters for external flow
  double MaInfinity, Aoa, PsInfinity, TsInfinity, RhoInfinity, velInfinity,
      uInfinity, vInfinity;

  // Physics parameters for internal flow
  double PtInlet, TtInlet, flowAngIn, PsOutlet, approxFlowAngOut, PsRatio;

  // Boundary Conditions
  std::unordered_map<std::string, BoundaryType> boundaryMap;
  std::vector<PeriodicInfo> periodics;
  std::set<std::string> periodicMaster;

  // Geometrical reference values
  double xRefPoint, yRefPoint, chord;

  // Iteration control
  size_t maxIteration, numberIterDump;
  double convTol;
  bool restart;

  // Numerical parameter
  ConvectionScheme convecScheme;
  double lusgsParameter;
  TemporalScheme temporalScheme;
  Limiter limiterType;
  double CFL, imResiSmooth;
  int numOfIterSmooth;
  timeStep timestep_;
  bool preCon;
  double preConFactor;
  int spatialOrder;
  double entropyCorreRoe, limiterCoeff;
  bool farCorrect;

  int temperalStages;

  std::vector<double> stageCoeff;
  std::vector<double> dissipationBlend;
  std::vector<int> dissipationEval;
};

}  // namespace preprocess