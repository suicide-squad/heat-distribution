cmake_minimum_required(VERSION 3.6)
project(RungeKutta)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c11")


set(SOURCE_FILES rugekutta.c)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})