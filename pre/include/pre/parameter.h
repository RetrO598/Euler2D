#pragma once
#include <cstddef>
#include <string>
#include <vector>

namespace preprocess {

enum class flowType { Internal, External };

enum class equationType { Euler, NavierStokes };

enum class timeStep { local, global };

struct parameter {
  explicit parameter(){};

  parameter(const parameter &other) = delete;

  parameter &operator=(const parameter &other) = delete;

  parameter(parameter &&other) = delete;

  parameter &operator=(parameter &&other) = delete;

  ~parameter() = default;

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

  // Geometrical reference values
  double xRefPoint, yRefPoint, chord;

  // Iteration control
  size_t maxIteration, numberIterDump;
  double convTol;
  bool restart;

  // Numerical parameter
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