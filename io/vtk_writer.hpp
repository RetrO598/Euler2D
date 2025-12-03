#pragma once

#include "mesh/mesh_data.hpp"
#include "solver/FVMSolver.h"

#include <string>

namespace io {

/**
 * @brief Writes the solution data to a VTK unstructured grid file (.vtu).
 *
 * This function generates a VTK file containing the mesh and cell-centered solution
 * data, including conservative variables, primitive variables (calculated), and pressure.
 *
 * @param filename The path to the output file (e.g., "solution.vtu").
 * @param mesh The computational mesh, containing points and element connectivity.
 * @param solver The FVM solver instance, containing the solution vectors.
 */
void writeVTKFile(const std::string& filename,
                  const preprocess::MeshData<2>& mesh,
                  const solver::FVMSolver& solver);

} // namespace io
