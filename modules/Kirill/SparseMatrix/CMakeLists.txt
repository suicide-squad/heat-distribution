cmake_minimum_required(VERSION 2.8)
project(SparseMatrix)

set(GCC_COVERAGE_COMPILE_FLAGS "-std=c++14")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

get_filename_component(DIR_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
set(MODULE      "${DIR_NAME}")
set(LIBRARY     "lib_${MODULE}")
set(TESTS       "test_${MODULE}")
set(APPLICATION "app_${MODULE}")

include_directories("./include")

add_subdirectory(src)
add_subdirectory(test)