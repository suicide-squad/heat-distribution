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
                valuesCounter++;

                if (pointerB.empty()) {
                    pointerB.push_back(j + i*j);
                } else if (pointerB.size() == i) {
                    pointerB.push_back(valuesCounter + pointerB[0] - 1);
                }
                columns.push_back(j);
            }
        }
        if (pointerE.empty()) {
            pointerE.push_back(valuesCounter + 1);
        } else if (pointerE.size() == i) {
            pointerE.push_back(valuesCounter + 1);
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
}