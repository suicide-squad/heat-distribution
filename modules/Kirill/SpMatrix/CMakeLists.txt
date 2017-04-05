cmake_minimum_required(VERSION 2.8)

#find_program(ICC_PATH  NAMES icc)
#find_program(ICPC_PATH NAMES icpc)
#message(${ICC_PATH})
#message(${ICPC_PATH})
#
#if(${ICC_PATH} MATCHES ICC_PATH-NOTFOUND OR
#        ${ICPC_PATH} MATCHES ICPC_PATH-NOTFOUND)
#    message(STATUS "Cannot find Intel compiler. You may need to run `. /opt/intel/bin/compilervars.sh intel64'")
#else()
#
#    set(CMAKE_C_COMPILER     ${ICC_PATH})
#    set(CMAKE_CXX_COMPILER ${ICPC_PATH})
#endif()

project(SpMatrix)

if (UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
endif(UNIX)

###################### OpenMP ###################################
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()
#################################################################

# math
find_library( M_LIB m )

######################### MPI ###################################
find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_PATH} )
#################################################################

add_definitions(${CMAKE_CXX_FLAGS})
add_definitions(${CMAKE_C_FLAGS})

get_filename_component(DIR_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
set(MODULE      "${DIR_NAME}")
set(LIBRARY     "lib_${MODULE}")
set(TESTS       "test_${MODULE}")
set(EULER       "euler_${MODULE}")
set(RUNGEKUTT   "runge_${MODULE}")
set(IMPLICIT    "implicit_euler_${MODULE}")

include_directories("./include")

add_subdirectory(src)
# add_subdirectory(test)
add_subdirectory(app)