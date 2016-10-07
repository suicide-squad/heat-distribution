cmake_minimum_required(VERSION 3.6)
project(HeatDistribution)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(MODULE  "HeatDistribution")
set(LIBRARY ${MODULE})
set(APPLICATION "app_${MODULE}")

include_directories("./include")

add_subdirectory(src)
add_subdirectory(app)
