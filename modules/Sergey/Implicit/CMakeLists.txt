set(target "${APP}_Implicit")

set(SOURCE_FILES main.cpp SparseMatrix.cpp SparseMatrix.h)
add_executable(${target} ${SOURCE_FILES})