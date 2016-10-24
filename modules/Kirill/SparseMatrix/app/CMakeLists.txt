set(target ${APPLICATION})

file(GLOB srcs "*.cpp")

add_executable(${target} ${srcs})

target_link_libraries(${target} ${LIBRARY})
