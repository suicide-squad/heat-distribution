cmake_minimum_required(VERSION 3.6)
project(SparseMatrix)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES_EULER src/main_Euler.cpp src/SparseMatrix.cpp src/SparseMatrix.h )
add_executable(${PROJECT_NAME}_Euler ${SOURCE_FILES_EULER})

set(SOURCE_FILES_RUNGE src/main_Runge.cpp src/SparseMatrix.cpp src/SparseMatrix.h )
add_executable(${PROJECT_NAME}_Runge ${SOURCE_FILES_RUNGE})