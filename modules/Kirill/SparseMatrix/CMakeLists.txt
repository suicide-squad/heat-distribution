cmake_minimum_required(VERSION 2.8)
project(SparseMatrix)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS}")
endif()

add_definitions(${CMAKE_CXX_FLAGS})

get_filename_component(DIR_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
set(MODULE      "${DIR_NAME}")
set(LIBRARY     "lib_${MODULE}")
set(TESTS       "test_${MODULE}")
set(EULER "euler_${MODULE}")
set(RUGEKUTT "runge_${MODULE}")

include_directories("./include")

add_subdirectory(src)
add_subdirectory(test)
#add_subdirectory(app)