set(target ${LIBRARY})

set( srcs methods.cpp euler.cpp )

add_library(${target} STATIC ${srcs})