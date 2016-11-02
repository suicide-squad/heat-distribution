//
// Created by lenferd on 27.10.16.
//

#include "SparseMatrix.h"


void SparseMatrix::fillMatrix(double** &matrix, int widthSize, int heightSize) {
    int valuesCounter = 0;
    for (int i = 0; i < heightSize; ++i ) { // line

        for (int j = 0; j < widthSize; ++j) {   // column
            if (matrix[i][j]!= 0 ) {
                values.push_back(matrix[i][j]);

                if (pointerB.empty()) {
                    pointerB.push_back(valuesCounter);
                } else if (pointerB.size() == i) {
                    pointerB.push_back(valuesCounter);
                }
                columns.push_back(j);
                valuesCounter++;
            }
        }
        if (pointerE.empty()) {
            pointerE.push_back(valuesCounter);
        } else if (pointerE.size() == i) {
            pointerE.push_back(valuesCounter);
        }
    }
}

void SparseMatrix::printVectors() {
    printf("values\n");
    for (auto x : values)
        printf("%lf ", x);
    printf("\n");

    printf("columns\n");
    for (auto x : columns)
        printf("%d ", x);
    printf("\n");

    printf("pointerB\n");
    for (auto x : pointerB)
        printf("%d ", x);
    printf("\n");

    printf("pointerE\n");
    for (auto x : pointerE)
        printf("%d ", x);
    printf("\n");
}

double* SparseMatrix::multiplicateVector(vector<double> vect) {
    double* resultVector = new double[vect.size()];
    for (int i = 0; i < vect.size(); i++) {
        resultVector[i] = 0;
    }

    int index = 0;
    for (int i = 0; i < vect.size(); i++){  // iteration FOR RESULT VECTOR!!!
        while (index < pointerE[i]) {
            resultVector[i] += values[index] * vect[columns[index]];
            ++index;
        }
    }

    return  resultVector;
}

void SparseMatrix::multiplicateVector(double *&vect, double *&result, int size) {
    int index = 0;

    for (int j = 0; j < size; ++j) {
        result[j] = 0;
    }

    for (int i = 0; i < size; i++){  // iteration FOR RESULT VECTOR!!!
        while (index < pointerE[i]) {
            result[i] += values[index] * vect[columns[index]];
            ++index;
        }
        result[i] += vect[i];
    }
}

void SparseMatrix::testEuler(int size, double expr) {
    int index = 0;
    values.push_back(1);
    columns.push_back(0);
    pointerB.push_back(index++);
    // TODO remove this array
    pointerE.push_back(index);
    for (int i = 1; i < size - 1; ++i) {

        values.push_back(expr);
        columns.push_back(i-1);
        pointerB.push_back(index++);

        values.push_back(-2*expr);
        columns.push_back(i);
        ++index;

        values.push_back(expr);
        columns.push_back(i+1);
        ++index;

        pointerE.push_back(index);
    }
    values.push_back(1);
    columns.push_back(size - 1);
    pointerB.push_back(index++);
    pointerE.push_back(index);
}