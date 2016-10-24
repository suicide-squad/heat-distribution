//
// Created by kirill on 23.10.16.
//

#ifndef SPARSE_SPARSE_H
#define SPARSE_SPARSE_H

#include <vector>

typedef double TYPE;
typedef std::vector<TYPE> vector;

// CSR (Compressed Sparse Rows)

class SpareMatrix {
 public:
  explicit SpareMatrix(TYPE*, int*, int*, const int, const int);
  ~SpareMatrix();
  vector operator*(const vector&);


 private:
  TYPE* value_;   // Элементы матрицы
  int* col_;      // Номера столбцов для каждого элемента
  int* rowIndex_; // Место каждого ненулевого элемента в каждой строке
  int nz_;        // Количество ненулевых
  int nRows_;     // Количество строк
};



#endif //SPARSE_SPARSE_H
