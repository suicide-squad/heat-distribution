set(target ${IMPLICIT})

file(GLOB SOURCE_DOUBLE_EULER "implicit_euler_sparse.c")
add_executable(${target} ${SOURCE_DOUBLE_EULER})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${M_LIB})