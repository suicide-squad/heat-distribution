//
// Created by lenferd on 27.10.16.
//

#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(int size, int rows) {
    _size = size;
    _rows = rows;
    values = new double[size];
    columns = new int[size];
    pointerB = new int[rows+1];
}


void SparseMatrix::multiplicateVector(double *&vect, double *&result, int size) {
    int index = 0;

    omp_set_num_threads(2);
    #pragma omp parallel for if (ENABLE_PARALLEL)
    for (int j = 0; j < size; ++j) {
        result[j] = 0;
    }

    #pragma omp parallel for if (ENABLE_PARALLEL)
    for (int i = 0; i < size; i++){  // iteration FOR RESULT VECTOR!!!
        while (index < pointerB[i+1]) {
            result[i] += values[index] * vect[columns[index]];
            ++index;
        }
    }
}

void SparseMatrix::fillMatrix2Expr(int size, double expr1, double expr2) {
    int index = 0;
    int pIndex = 0;

    values[index] = 1;
    columns[index] = 0;
    pointerB[pIndex++] = 0;
    ++index;

    for (int i = 1; i < size - 1; ++i) {

        //printf("index %d \n", index);
        values[index] = expr1;
        columns[index] = i - 1;
        pointerB[pIndex++] = index;
        ++index;

        values[index] = expr2;
        columns[index] = i;
        ++index;

        values[index] = expr1;
        columns[index] = i + 1;
        ++index;
    }

    values[index] = 1;
    columns[index] = size - 1;
    pointerB[pIndex++] = index;

    pointerB[pIndex] = index + 1;   //end
}

void SparseMatrix::printVectors() {
    printf("values\n");
    for (int i = 0; i < _size; ++i) {
        printf("%lf ", values[i]);
    }
    printf("\n");

    printf("columns\n");
    for (int i = 0; i < _size; ++i) {
        printf("%d ", columns[i]);
    }
    printf("\n");

    printf("pointerB\n");
    for (int i = 0; i < _rows + 1; ++i) {
        printf("%d ", pointerB[i]);
    }
    printf("\n");
}
