cmake_minimum_required(VERSION 3.6)
project(SparseMatrix)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES_EULER src/main_Euler.cpp src/SparseMatrix.cpp src/SparseMatrix.h )
add_executable(${PROJECT_NAME}_Euler ${SOURCE_FILES_EULER})

set(SOURCE_FILES_RUNGE src/main_Runge.cpp src/SparseMatrix.cpp src/SparseMatrix.h )
add_executable(${PROJECT_NAME}_Runge ${SOURCE_FILES_RUNGE})