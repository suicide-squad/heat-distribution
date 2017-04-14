#include_directories(include)
#
#get_filename_component(UNDER_DIR_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
#
#set(DIR_NAME "${DIR_NAME}_${UNDER_DIR_NAME}")
#set(MODULE      "${DIR_NAME}")
#set(APP         "app_${MODULE}")
#set(LIBRARY     "lib_${MODULE}")
#
#add_subdirectory(src)
#add_subdirectory(app)

set(target kekus)

find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_PATH} )
include_directories(include)

set(SOURCE_FILES app/Euler.cpp src/SparseMatrix.cpp include/SparseMatrix.h)
add_executable(${target} ${SOURCE_FILES})
target_link_libraries(${target} ${MPI_LIBRARIES})
