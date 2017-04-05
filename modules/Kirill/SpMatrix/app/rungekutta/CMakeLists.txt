set(target ${RUNGEKUTT})
set(target2 ${RUNGEKUTT}_complex)

file(GLOB SOURCE_DOUBLE_EULER "rungekutta_sparse.c")
add_executable(${target} ${SOURCE_DOUBLE_EULER})

target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_LIBRARIES})

# file(GLOB SOURCE_COMPLEX_EULER "rungekutta_sparse_complex.c")
# add_executable(${target2} ${SOURCE_COMPLEX_EULER})

# target_link_libraries(${target2} ${LIBRARY})