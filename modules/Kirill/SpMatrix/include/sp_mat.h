//
// Created by kirill on 23.10.16.
//

#ifndef SPARSE_SPARSE_H
#define SPARSE_SPARSE_H

#define ENABLE_PARALLEL 1
#define _COMPLEX_

#include <omp.h>

#ifdef __cplusplus
extern "C" {
#endif


#ifdef _COMPLEX_

#include <complex.h>
typedef complex TYPE;

#else
typedef double TYPE;
#endif

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


void sum(int N, double h, TYPE **result, TYPE *U, TYPE *k1, TYPE *k2, TYPE *k3, TYPE *k4);


#ifdef __cplusplus
}
#endif


#endif //SPARSE_SPARSE_H
