#include <stdio.h>
#include <stdlib.h>

#include "sparse_matrix.h"

void initSpMat(spMatrix* mat, int nz, int nRows) {
  mat->nz = nz;
  mat->nRows = nRows;
  mat->value = (TYPE *)malloc(sizeof(TYPE) * nz);
  mat->col = (int *)malloc(sizeof(int) * nz);
  mat->rowIndex = (int *)malloc(sizeof(int) * (nRows) + 1);
}

void freeSpMat(spMatrix* mat) {
  free(mat->value);
  free(mat->col);
  free(mat->rowIndex);
}

void multMV(TYPE** result, spMatrix mat, TYPE* vec) {
  TYPE localSum;
  #pragma omp parallel private(localSum) num_threads(2) if (ENABLE_PARALLEL)
  {
    #pragma omp for nowait
    for (int i = 0; i < mat.nRows; i++) {
      localSum = 0.0;
      for (int j = mat.rowIndex[i]; j < mat.rowIndex[i + 1]; j++)
        localSum += mat.value[j] * vec[mat.col[j]];
      (*result)[i] = localSum;
    }
  }
}

void sum(int N, double h, TYPE **result, TYPE *U, TYPE *k1, TYPE *k2, TYPE *k3, TYPE *k4) {
#pragma omp parallel for num_threads(2) if (ENABLE_PARALLEL)
  for (int i = 0; i < N; i++)
    (*result)[i] = U[i] + h*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}


//void SpareMatrix::print() {
//  for (int i = 0; i < nRows_; i++) {
//    for (int j = 0; j < nRows_; j++)
//      printf("%3.7lf\t", this->operator()(i, j));
//    printf("\n");
//  }
//}
//TYPE SpareMatrix::operator()(int i, int j) {
//  TYPE result = 0;
//  int N1 = rowIndex_[i];
//  int N2 = rowIndex_[i+1];
//  for(int k = N1; k < N2; k++)
//  {
//    if (col_[k] == j)
//    {
//      result = value_[k];
//      break;
//    }
//  }
//  return result;
//}
