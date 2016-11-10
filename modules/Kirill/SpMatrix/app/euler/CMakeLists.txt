set(target ${EULER})
set(target2 ${EULER}_complex)

file(GLOB SOURCE_DOUBLE_EULER "euler_sparse.c")
add_executable(${target} ${SOURCE_DOUBLE_EULER})
target_link_libraries(${target} ${LIBRARY})

file(GLOB SOURCE_COMPLEX_EULER "euler_sparse_comlex.c")
add_executable(${target2} ${SOURCE_COMPLEX_EULER})

target_link_libraries(${target2} ${LIBRARY})
