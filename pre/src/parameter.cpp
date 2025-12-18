#include <pre/parameter.h>

#include <iostream>

namespace preprocess {
void parameter::printParameters() {
  if (!initialized) [[unlikely]] {
    std::cout << "Parameters not initialized" << "\n";
    return;
  }
  std::cout << "---------------------------------" << "\n";
  std::cout << "Grid file: " << gridFile << "\n";
  std::cout << "Flow field to plot: " << flowFieldPlot << "\n";
  std::cout << "Surface quantities to plot: " << surfacePlot << "\n";
  std::cout << "Convergence history: " << convHistory << "\n";
  std::cout << "Restart solution (in): " << restartIn << "\n";
  std::cout << "Restart solution (out): " << restartOut << "\n";

  std::cout << "# Physics - general" << "\n";
  std::cout << "------------------------------------" << "\n";
  std::cout << "Flow Type: ";
  if (flowtype_ == flowType::External) {
    std::cout << "External flow" << "\n";
  } else if (flowtype_ == flowType::Internal) {
    std::cout << "Internal flow" << "\n";
  }

  std::cout << "Equation Type: ";
  if (equationtype_ == equationType::Euler) {
    std::cout << "Euler equation" << "\n";
  } else if (equationtype_ == equationType::NavierStokes) {
    std::cout << "Navier-Stokes equation" << "\n";
  } else if (equationtype_ == equationType::RANS) {
    std::cout << "RANS equation" << "\n";
  }

  std::cout << "Ratio of specific heats: " << gamma << "\n";
  std::cout << "Specific heat coeff. at constant pressure: " << Cp << "\n";
  std::cout << "Reynolds number: " << Re << "\n";
  std::cout << "Reference velocity: " << refVel << "\n";
  std::cout << "Reference Density: " << refRho << "\n";
  std::cout << "Laminar viscosity: " << refVisc << "\n";
  std::cout << "Laminar Prandtl number: " << Prandtl << "\n";

  if (flowtype_ == flowType::External) {
    std::cout << "# Physics - external flow" << "\n";
    std::cout << "--------------------------------------" << "\n";
    std::cout << "Mach-number at infinity: " << MaInfinity << "\n";
    std::cout << "Angle of attack: " << Aoa << "\n";
    std::cout << "Static pressure at infinity: " << PsInfinity << "\n";
    std::cout << "Static temperature at infinity: " << TsInfinity << "\n";
  } else {
    std::cout << "# Physics - internal flow" << "\n";
    std::cout << "--------------------------------------" << "\n";
    std::cout << "Total pressure at inlet: " << PtInlet << "\n";
    std::cout << "Total temperature at inlet: " << TtInlet << "\n";
    std::cout << "Flow angle at inlet (with x-axis): " << flowAngIn << "\n";
    std::cout << "Static pressure at outlet: " << PsOutlet << "\n";
    std::cout << "Approximate flow angle at outlet (with x-axis): "
              << approxFlowAngOut << "\n";
    std::cout << "Approximate ratio of inlet to outlet static pressure: "
              << PsRatio << "\n";
  }

  std::cout << "# Geometrical reference values" << "\n";
  std::cout << "----------------------------------------" << "\n";
  std::cout << "x-coordinate of reference point: " << xRefPoint << "\n";
  std::cout << "y-coordinate of reference point: " << yRefPoint << "\n";
  std::cout << "Reference or chord length: " << chord << "\n";

  std::cout << "# Iteration control" << "\n";
  std::cout << "----------------------------------------" << "\n";
  std::cout << "Max number of iteration: " << maxIteration << "\n";
  std::cout << "Number of iterations between solution dumps: " << numberIterDump
            << "\n";
  std::cout << "Convergence tolerance: " << convTol << "\n";
  std::cout << "Use previous solution for restart: " << restart << "\n";

  std::cout << "# Numerical parameters" << "\n";
  std::cout << "----------------------------------------" << "\n";
  std::cout << "CFL number: " << CFL << "\n";
  std::cout << "Coefficient of implicit residual smoothing: " << imResiSmooth
            << "\n";
  std::cout << "Number of Jacobi iterations for residual smoothing: "
            << numOfIterSmooth << "\n";
  std::cout << "Local/Global time step: ";
  if (timestep_ == timeStep::local) {
    std::cout << "local" << "\n";
  } else if (timestep_ == timeStep::global) {
    std::cout << "global" << "\n";
  }
  std::cout << "Low Mach number preconditioning: " << preCon << "\n";
  std::cout << "Precondition parameter: " << preConFactor << "\n";
  std::cout << "First-order / Second-order: " << spatialOrder << "\n";
  std::cout << "Limiter coefficient (used only for 2nd-order Roe scheme): "
            << limiterCoeff << "\n";
  std::cout << "Entropy correction coefficient: " << entropyCorreRoe << "\n";
  std::cout << "Correction of far-field due to single vortex: " << farCorrect
            << "\n";
  std::cout << "Number of stages: " << temperalStages << "\n";
  std::cout << "Stage coefficients: ";
  for (const auto& i : stageCoeff) {
    std::cout << i << " ";
  }
  std::cout << "\n";
  std::cout << "Dissipation blending coeff.: ";
  for (const auto& i : dissipationBlend) {
    std::cout << i << " ";
  }
  std::cout << "\n";

  std::cout << "Dissipation  evaluation (1=yes): ";
  for (const auto& i : dissipationEval) {
    std::cout << i << " ";
  }
  std::cout << "\n";

  if (flowtype_ == flowType::External) {
    std::cout << "Density at infinity: " << RhoInfinity << "\n";
    std::cout << "Velocity at infinity: " << velInfinity << "\n";
    std::cout << "u at infinity: " << uInfinity << "\n";
    std::cout << "v at infinity: " << vInfinity << "\n";
    std::cout << "Reference Mach pow: " << refMach2 << "\n";
    std::cout << "Reference Density: " << refRho << "\n";
    std::cout << "Reference Velocity: " << refVel << "\n";
  } else {
    std::cout << "Reference Mach pow: " << refMach2 << "\n";
  }
}
}  // namespace preprocess