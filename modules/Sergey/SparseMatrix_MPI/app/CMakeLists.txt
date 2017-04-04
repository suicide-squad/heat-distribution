# Эйлер
set(target "${APP}_Euler")

file(GLOB src_euler Euler.cpp)

add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})

#set(target "${APP}_SparseMatrix_MPI_Euler")
#
#set(SOURCE_FILES_EULER
#        src/SparseMatrix.cpp
#        )
#
#add_executable(${target} ${SOURCE_FILES_EULER})
#target_link_libraries(${target} ${MPI_LIBRARIES})