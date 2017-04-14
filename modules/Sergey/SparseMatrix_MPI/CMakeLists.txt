cmake_minimum_required(VERSION 3.6)
#project(SparseMatrixMPI)
set(target Sergey_SparseMatrix_MPI_Euler)

find_package(MPI REQUIRED)
include_directories( ${MPI_INCLUDE_PATH} )

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${OpenMP_C_FLAGS}")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}    ${OpenMP_C_FLAGS}")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp SparseMatrix.cpp SparseMatrix.h)
add_executable(${target} ${SOURCE_FILES})
target_link_libraries(${target} ${MPI_LIBRARIES})
