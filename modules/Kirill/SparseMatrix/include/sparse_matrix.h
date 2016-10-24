//
// Created by kirill on 23.10.16.
//

#ifndef SPARSE_SPARSE_H
#define SPARSE_SPARSE_H

#include <omp.h>

#include <vector>
#include <algorithm>

typedef double TYPE;
typedef std::vector<TYPE> vector;

#define ENABLE_PARALLEL 1

// CSR (Compressed Sparse Rows)
class SpareMatrix {
 public:
  explicit SpareMatrix(TYPE*, int*, int*, const size_t , const size_t);
  ~SpareMatrix();
  vector operator*(const vector&);


 private:
  TYPE* value_;   // Элементы матрицы
  int* col_;      // Номера столбцов для каждого элемента
  int* rowIndex_; // Место каждого ненулевого элемента в каждой строке
  size_t nz_;        // Количество ненулевых
  size_t nRows_;     // Количество строк
};



#endif //SPARSE_SPARSE_H
