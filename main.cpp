#include <pre/geometry.h>
#include <pre/parameter.h>
#include <pre/reader.h>
#include <solver/FVMSolver.h>

#include <chrono>
#include <grid/grid_builder.hpp>
#include <iostream>
#include <string>

#include "io/vtk_writer.hpp"
#include "mesh/mesh_data.hpp"
#include "mesh/mesh_reader.hpp"
#include "solver/timeIntegrator.hpp"

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;
  }
  std::string inputFile = argv[1];
  preprocess::yamlReader reader(inputFile);
  preprocess::parameter param;
  // reader.customRead(param);
  reader.read(param);
  param.printParameters();

  preprocess::Geometry geometry(param);
  std::unique_ptr<preprocess::MeshDataBase> mesh;

  std::string extension = ".su2";
  if (param.gridFile.size() >= extension.size() &&
      param.gridFile.compare(param.gridFile.size() - extension.size(),
                             extension.size(), extension) == 0) {
    // geometry.ReadSU2Grid();
    // geometry.printInfo();
    // geometry.ComputeMetrics();
    std::cout << "reading SU2 format mesh" << "\n";
    preprocess::SU2Reader reader(param.gridFile);
    mesh = reader.readMesh();
    auto &data = static_cast<preprocess::MeshData<2> &>(*mesh);
    preprocess::GridBuilder::FVMBuilder::build(data);
    geometry = data.changeTo2Dgeometry(param.boundaryMap);
    geometry.printInfo();
  } else {
    geometry.ReadGrid();
    geometry.printInfo();
    geometry.ComputeMetrics();
  }

  solver::FVMSolver solver(param, geometry);
  solver.initSolver();
  solver.ConvToDependAll();
  solver.BoundaryConditions();

  solver.limiter->limiterRefVals();
  solver.iter = 0;

  if (param.temporalScheme == preprocess::TemporalScheme::LUSGS) {
    std::cout << "using LUSGS temporal scheme." << "\n";
    solver::LUSGSIntegrator integrator(solver);
    solver.computeWaveSpeed();
    auto start = std::chrono::high_resolution_clock::now();
    do {
      solver.iter++;
      integrator.timeAdvance();
      solver.Convergence();
    } while (!solver.Converged());
    auto end = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration<double>(end - start);
    std::cout << "Computation costs: " << time.count() << "s" << "\n";
  } else if (param.temporalScheme == preprocess::TemporalScheme::RungeKutta) {
    std::cout << "using RungeKutta temporal scheme." << "\n";
    solver::RungeKuttaTimeIntegrator integrator(solver);

    auto start = std::chrono::high_resolution_clock::now();
    do {
      solver.iter++;
      integrator.timeAdvance();
      solver.Convergence();
    } while (!solver.Converged());
    auto end = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration<double>(end - start);
    std::cout << "Computation costs: " << time.count() << "s" << "\n";
  }

  // solver.computeWaveSpeed();
  // do {
  //   solver.iter++;
  //   integrator.timeAdvance();
  //   // solver.solve();
  //   solver.Convergence();
  // } while (!solver.Converged());

  // solver.writeTecplotDat();
  // solver.writeLineDat();

  if (mesh) {
    std::cout << "Writing VTK output file..." << std::endl;
    auto &meshData = static_cast<preprocess::MeshData<2> &>(*mesh);
    std::string vtkFilename = param.title + ".vtu";
    io::writeVTKFile(vtkFilename, meshData, solver);
  }

  return 0;
}