#include <stdio.h>
#include <stdlib.h>
#include <printf.h>

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


//SpareMatrix::SpareMatrix(TYPE* value, int* col, int* rowIndex, const size_t nz, const size_t nRows) {
//  nz_ = nz;
//  nRows_ = nRows;
//  value_ = new TYPE[nz_];
//  col_ = new int[nz_];
//  rowIndex_ = new int[nRows_ + 1];
//  std::copy(value, value + nz, value_);
//  std::copy(col, col + nz, col_);
//  std::copy(rowIndex, rowIndex + nRows + 1, rowIndex_);
//}
//
//SpareMatrix::~SpareMatrix() {
//  delete[] value_;
//  delete[] col_;
//  delete[] rowIndex_;
//}
//



//vector SpareMatrix::operator*(const vector& v) {
//  vector result(nRows_);
//  TYPE localSum;
//  #pragma omp parallel private(localSum) num_threads(2) if (ENABLE_PARALLEL)
//  {
//    #pragma omp for nowait
//    for (int i = 0; i < nRows_; i++) {
//      localSum = 0;
//      for (int j = rowIndex_[i]; j < rowIndex_[i + 1]; j++)
//        localSum += value_[j] * v[col_[j]];
//      result[i] = localSum;
//    }
//  }
//  return result;
//}
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
