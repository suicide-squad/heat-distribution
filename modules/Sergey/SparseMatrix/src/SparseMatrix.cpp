//
// Created by lenferd on 27.10.16.
//

#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(double** &matrix, int widthSize, int heightSize) {
    int valuesCounter = 0;
    for (int i = 0; i < heightSize; ++i ) {
        for (int j = 0; j < widthSize; ++j) {
            if (matrix[i][j]!= 0 ) {
                values.push_back(matrix[i][j]);
                valuesCounter++;
                if (colums.empty()) {
                    pointerB.push_back(valuesCounter);
                } else {
                    pointerB.push_back(valuesCounter + pointerB[0] - 1);
                }
                colums.push_back(i);
            }

        }
        if (i == 1) {
            pointerE.push_back(valuesCounter);
        /*} else if (i = heightSize - 1 ) {
            pointerE.push_back(valuesCounter);*/
        } else {
            pointerE.push_back(valuesCounter + pointerE[0] - 1);
        }
    }
}