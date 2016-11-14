//
// Created by lenferd on 27.10.16.
//

#ifndef SPARSEMATRIX_SPARSEMATRIX_H
#define SPARSEMATRIX_SPARSEMATRIX_H
#include <omp.h>
#include <cstdio>


const int ENABLE_PARALLEL = 0;


class SparseMatrix {
private:
    int _size;
    int _rows;
    double* values;
    int* columns;
    int* pointerB;

public:
    SparseMatrix(int size, int rows) ;
    void multiplicateVector(double* &vect, double* &result, int size);
    void fillMatrix2Expr(int size, double expr1, double expr2);
    void printVectors();
};


#endif //SPARSEMATRIX_SPARSEMATRIX_H
