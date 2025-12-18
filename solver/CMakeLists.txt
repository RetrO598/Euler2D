file(GLOB_RECURSE srcs CONFIGURE_DEPEnDS src/*.cpp include/*.h include/*.hpp)
add_library(solver STATIC ${srcs})
target_include_directories(solver PUBLIC include)

target_link_libraries(solver PUBLIC preprocess)
