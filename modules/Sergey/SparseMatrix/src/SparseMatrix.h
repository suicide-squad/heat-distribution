//
// Created by lenferd on 27.10.16.
//

#ifndef SPARSEMATRIX_SPARSEMATRIX_H
#define SPARSEMATRIX_SPARSEMATRIX_H
#include <vector>

using std::vector;

class SparseMatrix {
    vector<double> values = new vector<double >;
    vector<double> colums;
    vector<double> pointerB;
    vector<double> pointerE;
    /*
    double* values = nullptr;
    double* columns = nullptr;
    double* pointerB = nullptr;
    double* pointerE = nullptr;
    */
    SparseMatrix(double** &matrix, int widthSize, int heightSize);
    SparseMatrix(double* &matrix, int size);
};


#endif //SPARSEMATRIX_SPARSEMATRIX_H
