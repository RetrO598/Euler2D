# Find the VTK package and its required components
find_package(VTK COMPONENTS 
  CommonCore
  IOXML
REQUIRED)

# Include VTK's headers
include(${VTK_USE_FILE})

# Define a library for our IO module
add_library(io vtk_writer.cpp)

# Link our library against the necessary VTK modules
target_link_libraries(io PRIVATE 
  VTK::CommonCore
  VTK::IOXML
)

# Add include directories to allow our code to find the VTK, pre, and solver headers
target_include_directories(io PUBLIC 
  "${PROJECT_SOURCE_DIR}/pre/include"
  "${PROJECT_SOURCE_DIR}/solver/include"
)
