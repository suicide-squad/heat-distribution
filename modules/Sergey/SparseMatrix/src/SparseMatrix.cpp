//
// Created by lenferd on 27.10.16.
//

#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(double **&matrix, int widthSize, int heightSize) {
    for (int i = 0; i < heightSize; ++i ) {
        for (int j = 0; j < widthSize; ++j) {
            if (matrix[i][j]!= 0 ) {
                values.push_back(matrix[i][j]);
                colums.push_back(i);
            }

        }
    }
}