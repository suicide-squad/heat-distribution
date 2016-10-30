//
// Created by kirill on 23.10.16.
//

#ifndef SPARSE_SPARSE_H
#define SPARSE_SPARSE_H

#include <omp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double TYPE;

#define ENABLE_PARALLEL 1

// CSR (Compressed Sparse Rows)
typedef struct {
  TYPE* value;   // Элементы матрицы
  int* col;      // Номера столбцов для каждого элемента
  int* rowIndex; // Место каждого ненулевого элемента в каждой строке
  int nz;        // Количество ненулевых
  int nRows;     // Количество строк
} spMatrix;

void initSpMat(spMatrix* mat, int nz, int nRows);
void freeSpMat(spMatrix* mat);

void multMV(TYPE** result, spMatrix matrix, TYPE* vector);

#ifdef __cplusplus
}
#endif


#endif //SPARSE_SPARSE_H
