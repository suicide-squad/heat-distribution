set(target ${LIBRARY})

set( srcs methods.cpp euler.cpp rungekutta.cpp )

add_library(${target} STATIC ${srcs})