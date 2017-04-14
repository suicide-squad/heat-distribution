# Эйлер
set(target "${APP}_Euler")
find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_PATH} )

file(GLOB src_euler Euler.cpp)

add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_LIBRARIES})
