set(target ${LIBRARY})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c11")

add_definitions(${CMAKE_CXX_FLAGS})

file(GLOB hrcs "${CMAKE_HOME_DIRECTORY}/include/*.h")
file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

add_library( ${target} STATIC ${hrcs} ${srcs} )
