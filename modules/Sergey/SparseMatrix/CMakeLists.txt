set(target "${APP}_SparseMatrix_Euler")

set(SOURCE_FILES_EULER src/main_Euler.cpp src/SparseMatrix.cpp src/SparseMatrix.h )
add_executable(${target} ${SOURCE_FILES_EULER})


set(target "${APP}_SparseMatrix_Runge")

set(SOURCE_FILES_RUNGE src/main_Runge.cpp src/SparseMatrix.cpp src/SparseMatrix.h )
add_executable(${target} ${SOURCE_FILES_RUNGE})