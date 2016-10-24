//
// Created by kirill on 23.10.16.
//

#ifndef SPARSE_SPARSE_H
#define SPARSE_SPARSE_H


// CSR (Compressed Sparse Rows)

typedef double TYPE;

class SpareMatrix {
 public:
  explicit SpareMatrix(TYPE*, int*, int*, const int, const int);
  ~SpareMatrix();
  TYPE*& operator*(TYPE *);


 private:
  TYPE* value_;   // Элементы матрицы
  int* col_;      // Номера столбцов для каждого элемента
  int* rowIndex_; // Место каждого ненулевого элемента в каждой строке
  int nz_;        // Количество ненулевых
  int nRows_;     // Количество строк
};



#endif //SPARSE_SPARSE_H
