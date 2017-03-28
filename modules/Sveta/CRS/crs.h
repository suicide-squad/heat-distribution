#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define ENABLE_PARALLEL 1
typedef struct {
    double* Value; //элементы матрицы
    int* col;//номера столбцов для ненулевых элементов
    int* rowindex;//индекс строки
    int nz;//количество ненулевых элементов
    int nRows;
} CRSMatrix;

void initCRSMartix(int nRows, int nz, CRSMatrix *crsm);
void freeCRSMatrix(CRSMatrix *);
void multCRSMatrix(double** result, CRSMatrix crsm, double* vec);
