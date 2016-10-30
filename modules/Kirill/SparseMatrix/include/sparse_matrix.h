//
// Created by kirill on 23.10.16.
//

#ifndef SPARSE_SPARSE_H
#define SPARSE_SPARSE_H

#include <omp.h>

#include <vector>
#include <algorithm>
#include <iostream>

typedef double TYPE;

class vector : public std::vector<TYPE> {
 public:
  vector(size_t size) : std::vector<TYPE>(size) {}
  vector() : std::vector<TYPE>() {}
  vector(std::initializer_list<TYPE> list) : std::vector<TYPE>(list) {}
  vector operator+ (const vector& b) {
    vector result(size());
    for(int i = 0; i < size(); ++i)
      result[i] = this->at(i) + b[i];
    return result;
  }
  void print() const {
    for (auto i = begin(); i!= end(); ++i)
      printf("%2.2lf \n", *i);
    printf("\n");
  }


  vector operator* (const TYPE b) {
    vector result(size());
    for (int i = 0; i < size(); ++i)
      result[i] = this->at(i)*b;
    return result;
  }
};

#define ENABLE_PARALLEL 1

// CSR (Compressed Sparse Rows)
class SpareMatrix {
 public:
  explicit SpareMatrix(TYPE*, int*, int*, const size_t , const size_t);
  ~SpareMatrix();
  vector operator*(const vector&);
  void print();
  TYPE operator()(int i, int j);



 private:
  TYPE* value_;   // Элементы матрицы
  int* col_;      // Номера столбцов для каждого элемента
  int* rowIndex_; // Место каждого ненулевого элемента в каждой строке
  size_t nz_;        // Количество ненулевых
  size_t nRows_;     // Количество строк
};



#endif //SPARSE_SPARSE_H
