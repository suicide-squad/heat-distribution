cmake_minimum_required(VERSION 3.6)
project(RungeKutta)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

add_definitions(${CMAKE_C_FLAGS})

set(SOURCE_FILES rungekutta.c)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})