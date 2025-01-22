file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h)
add_library(pre STATIC ${srcs})
target_include_directories(pre PUBLIC include)
