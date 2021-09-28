cmake_minimum_required(VERSION 3.12)
project(BlockEpsilonTree)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif ()

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    link_libraries(OpenMP::OpenMP_CXX)
endif ()

include_directories(lib/sdsl-lite/include include)

add_executable(example example.cpp)
