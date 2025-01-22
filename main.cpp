#include <pre/geometry.h>
#include <pre/parameter.h>
#include <pre/reader.h>
#include <solver/FVMSolver.h>

#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
  }
  std::string inputFile = argv[1];
  preprocess::simpleReader reader(inputFile);
  preprocess::parameter param;
  reader.customRead(param);
  param.printParameters();

  reader.close();

  preprocess::Geometry geometry(param.gridFile);

  geometry.ReadGrid();
  geometry.printInfo();
  geometry.ComputeMetrics();
  geometry.outputMeshInfo();
  solver::FVMSolver solver(param, geometry);
  solver.initSolver();
  solver.ConvToDependAll();
  solver.BoundaryConditions();

  solver.limiter->limiterRefVals(solver.param, solver.geom);
  solver.iter = 0;
  do {
    solver.iter++;
    solver.solve();
    solver.Convergence();
  } while (!solver.Converged());

  solver.writeTecplotDat();

  return 0;
}