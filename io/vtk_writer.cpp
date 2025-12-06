#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <iostream>
#include <vector>

#include "grid/point.hpp"
#include "pre/parameter.h"
#include "vtk_writer.hpp"

// Forward declaration needed for the friend declaration in FVMSolver.h
namespace solver {
class FVMSolver;
}

namespace io {

// This function needs access to private members of FVMSolver (cv and
// timeSteps). A 'friend' declaration must be added to the FVMSolver class
// definition.
void writeVTKFile(const std::string& filename,
                  const preprocess::MeshData<2>& mesh,
                  const solver::FVMSolver& solver) {
  auto points = vtkSmartPointer<vtkPoints>::New();
  auto cells = vtkSmartPointer<vtkCellArray>::New();
  auto cellTypes = vtkSmartPointer<vtkUnsignedCharArray>::New();
  auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // 1. Populate vtkPoints from the mesh's PointList.
  for (const auto& pt : mesh.PointList) {
    const auto& coords = pt.getCoord();
    points->InsertNextPoint(coords[0], coords[1], 0.0);
  }

  // 2. Populate vtkCellArray and cell types array for a mixed grid.
  for (const auto& elemHandle : mesh.ElemList) {
    const int vtkType = elemHandle.getVTKtype();
    const int numNodes = elemHandle.getnNodes();

    cellTypes->InsertNextValue(vtkType);

    auto vtk_ids = vtkSmartPointer<vtkIdList>::New();
    vtk_ids->SetNumberOfIds(numNodes);
    for (int i = 0; i < numNodes; ++i) {
      vtk_ids->SetId(i, elemHandle.getNodes(i));
    }
    cells->InsertNextCell(vtk_ids);
  }

  // 3. Assemble the unstructured grid.
  unstructuredGrid->SetPoints(points);
  unstructuredGrid->SetCells(cellTypes, cells);

  // 4. Prepare node-centered solution data.
  vtkPointData* pointData = unstructuredGrid->GetPointData();
  const size_t numPoints = mesh.PointList.size();

  // Access private members of FVMSolver (requires friendship).
  const std::vector<solver::CONS_VAR>& consVars = solver.cv;
  const std::vector<double>& timeSteps = solver.timeSteps;

  if (consVars.size() != numPoints || timeSteps.size() != numPoints) {
    std::cerr << "VTK Writer Error: Mismatch between number of points ("
              << numPoints << ") and solution data size (" << consVars.size()
              << ")." << std::endl;
    return;
  }

  // Create VTK data arrays
  auto density = vtkSmartPointer<vtkDoubleArray>::New();
  density->SetName("Density");
  density->SetNumberOfTuples(numPoints);

  auto momentum = vtkSmartPointer<vtkDoubleArray>::New();
  momentum->SetName("Momentum");
  momentum->SetNumberOfComponents(3);
  momentum->SetNumberOfTuples(numPoints);

  auto energy = vtkSmartPointer<vtkDoubleArray>::New();
  energy->SetName("Energy");
  energy->SetNumberOfTuples(numPoints);

  auto velocity = vtkSmartPointer<vtkDoubleArray>::New();
  velocity->SetName("Velocity");
  velocity->SetNumberOfComponents(3);
  velocity->SetNumberOfTuples(numPoints);

  auto pressure = vtkSmartPointer<vtkDoubleArray>::New();
  pressure->SetName("Pressure");
  pressure->SetNumberOfTuples(numPoints);

  auto timeStepArr = vtkSmartPointer<vtkDoubleArray>::New();
  timeStepArr->SetName("TimeStep");
  timeStepArr->SetNumberOfTuples(numPoints);

  auto mueArr = vtkSmartPointer<vtkDoubleArray>::New();
  if (solver.param.equationtype_ == preprocess::equationType::NavierStokes) {
    mueArr->SetName("dynamicViscosity");
    mueArr->SetNumberOfTuples(numPoints);
  }

  // Populate the arrays with data from the solver
  for (size_t i = 0; i < numPoints; ++i) {
    const auto& cv = consVars[i];
    const double rho = cv.dens;
    const double u = (rho == 0) ? 0 : cv.xmom / rho;
    const double v = (rho == 0) ? 0 : cv.ymom / rho;
    const double p = (solver.param.gamma - 1.0) *
                     (cv.ener - 0.5 * (cv.xmom * u + cv.ymom * v));

    density->SetTuple1(i, rho);
    momentum->SetTuple3(i, cv.xmom, cv.ymom, 0.0);
    energy->SetTuple1(i, cv.ener);
    velocity->SetTuple3(i, u, v, 0.0);
    pressure->SetTuple1(i, p);
    timeStepArr->SetTuple1(i, timeSteps[i]);

    if (solver.param.equationtype_ == preprocess::equationType::NavierStokes) {
      mueArr->SetTuple1(i, solver.dvlam[i].mu);
    }
  }

  // Add arrays to the grid's point data
  pointData->AddArray(density);
  pointData->AddArray(momentum);
  pointData->AddArray(energy);
  pointData->AddArray(velocity);
  pointData->AddArray(pressure);
  pointData->AddArray(timeStepArr);
  if (solver.param.equationtype_ == preprocess::equationType::NavierStokes) {
    pointData->AddArray(mueArr);
  }

  // 5. Write the grid to a .vtu file
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(unstructuredGrid);
  writer->SetDataModeToAscii();  // Or SetDataModeToBinary() for smaller files
  writer->Write();

  std::cout << "VTK file written: " << filename << std::endl;
}

}  // namespace io
