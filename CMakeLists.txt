cmake_minimum_required(VERSION 3.11)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXTENSIONS OFF)

project(CppCFD LANGUAGES CXX)

add_subdirectory(preprocess)
add_subdirectory(solver)
add_subdirectory(io)

add_executable(main main.cpp)
target_link_libraries(main PUBLIC preprocess solver io)
