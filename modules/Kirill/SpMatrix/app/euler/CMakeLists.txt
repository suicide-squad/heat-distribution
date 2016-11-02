set(target ${EULER})

file(GLOB srcs "*.c")

add_executable(${target} ${srcs})

target_link_libraries(${target} ${LIBRARY})
