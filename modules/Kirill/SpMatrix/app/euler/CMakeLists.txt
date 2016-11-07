set(target ${EULER})

file(GLOB SOURCE_DOUBLE_EULER "euler_sparse.c")
add_executable(${target} ${SOURCE_DOUBLE_EULER})

file(GLOB SOURCE_COMPLEX_EULER "euler_complex_sparse.c")
add_executable(${target}_complex ${SOURCE_COMPLEX_EULER})

target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target}_complex ${LIBRARY})
