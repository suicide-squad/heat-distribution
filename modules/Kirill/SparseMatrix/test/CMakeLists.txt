set(target ${TESTS})

file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable( ${target} ${srcs} )

target_link_libraries(${target} ${GTEST_BOTH_LIBRARIES} pthread)
target_link_libraries(${target} ${LIBRARY})