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
  float gamma, Cp, Re, refVel, refRho, Prandtl, refVisc, refMach2;

  // Physics parameters for external flow
  float MaInfinity, Aoa, PsInfinity, TsInfinity, RhoInfinity, velInfinity,
      uInfinity, vInfinity;

  // Physics parameters for internal flow
  float PtInlet, TtInlet, flowAngIn, PsOutlet, approxFlowAngOut, PsRatio;

  // Geometrical reference values
  float xRefPoint, yRefPoint, chord;

  // Iteration control
  size_t maxIteration, numberIterDump;
  float convTol;
  bool restart;

  // Numerical parameter
  float CFL, imResiSmooth;
  int numOfIterSmooth;
  timeStep timestep_;
  bool preCon;
  float preConFactor;
  int spatialOrder;
  float entropyCorreRoe, limiterCoeff;
  bool farCorrect;

  int temperalStages;

  std::vector<float> stageCoeff;
  std::vector<float> dissipationBlend;
  std::vector<int> dissipationEval;
};

}  // namespace preprocess