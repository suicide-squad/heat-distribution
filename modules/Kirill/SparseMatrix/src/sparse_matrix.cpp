#include "sparse_matrix.h"

#include <algorithm>

SpareMatrix::SpareMatrix(TYPE* value, int* col, int* rowIndex, const int nz, const int nRows) {
  nz_ = nz;
  nRows_ = nRows;
  value_ = new TYPE[nz_];
  col_ = new int[nz_];
  rowIndex_ = new int[nRows_ + 1];
  std::copy(value, value + nz, value_);
  std::copy(col, col + nz, col_);
  std::copy(rowIndex, rowIndex + nRows + 1, rowIndex_);
}

SpareMatrix::~SpareMatrix() {
  delete[] value_;
  delete[] col_;
  delete[] rowIndex_;
}

