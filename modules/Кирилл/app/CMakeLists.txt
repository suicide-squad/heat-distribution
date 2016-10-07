set(target ${APPLICATION})

set( SOURCE_FILES main.cpp )

add_executable(${target} ${SOURCE_FILES})
target_link_libraries(${target} ${LIBRARY})