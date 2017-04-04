//
// Created by lenferd on 27.10.16.
//

#include "SparseMatrix.h"


void spMatrixInit(SparseMatrix &sp, int size, int rows) {
    sp._size = size;
    sp._rows = rows;
    sp.values = new double[size];
    sp.columns = new int[size];
    sp.pointerB = new int[rows+1];
}

void multiplicateVector(SparseMatrix &sp, double *&vect, double *&result, int size) {

    omp_set_num_threads(4);

    #pragma omp parallel for if (ENABLE_PARALLEL)
    for (int i = 0; i < size; i++){  // iteration FOR RESULT VECTOR!!!
        double local_result = 0;
        for (int j = sp.pointerB[i]; j < sp.pointerB[i+1]; j++) {
            local_result += sp.values[j] * vect[sp.columns[j]];
        }
        result[i] = local_result;
    }
}

void fillMatrix2Expr(SparseMatrix &sp, int size, double expr1, double expr2) {
    int index = 0;
    int pIndex = 0;

    /** Boundaries rule
     *  If we on the edge, we should use same expression (line with parametrs), as the line after.
     *  If it's first line the pattern for her is line two. (+1)
     *  If it's last line, pattern - previous line.         (-1)
     *  Realization - fixes value, whose start to work, if we on the boundaries, joins @var x
     *  @var fixBounds
     */

    int fixBounds = 0;

    for (int i = 0; i < size; ++i) {
        if (i == 0 ) {
            fixBounds = 1;
        } else if ((i + 1) == size) {
            fixBounds = -1;
        }
        //printf("index %d \n", index);
        sp.values[index] = expr1;
        sp.columns[index] = fixBounds + i - 1;
        sp.pointerB[pIndex++] = index;
        ++index;

        sp.values[index] = expr2;
        sp.columns[index] = fixBounds + i;
        ++index;

        sp.values[index] = expr1;
        sp.columns[index] = fixBounds + i + 1;
        ++index;

        fixBounds = 0;
    }

    sp.pointerB[pIndex] = index + 1;   //end
}


void printVectors(SparseMatrix &sp) {
    printf("values\n");
    for (int i = 0; i < sp._size; ++i) {
        printf("%lf ", sp.values[i]);
    }
    printf("\n");

    printf("columns\n");
    for (int i = 0; i < sp._size; ++i) {
        printf("%d ", sp.columns[i]);
    }
    printf("\n");

    printf("pointerB\n");
    for (int i = 0; i < sp._rows + 1; ++i) {
        printf("%d ", sp.pointerB[i]);
    }
    printf("\n");
}
