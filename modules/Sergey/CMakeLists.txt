get_filename_component(DIR_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

set(MODULE      "${DIR_NAME}")
set(APP         "app_${MODULE}")

add_subdirectory(Euler)
add_subdirectory(Implicit)
add_subdirectory(Runge–Kutt)
add_subdirectory(SparseMatrix)
add_subdirectory(SparseMatrix_MPI)
