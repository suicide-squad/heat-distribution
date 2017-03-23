#include "crs.h"

void initCRSMartix(int nRows, int nz, CRSMatrix *crsm){
    crsm->nRows = nRows;
    crsm->nz = nz;
    crsm->Value = (double *)malloc(sizeof(double) * nz);
    crsm->col = (int *)malloc(sizeof(int) * nz);
    crsm->rowindex = (int *)malloc(sizeof(int) * (nRows + 1));
};
void freeCRSMatrix(CRSMatrix *crsm){
    free(crsm->Value);
    free(crsm->col);
    free(crsm->rowindex);
};
void multCRSMatrix(double** result, CRSMatrix crsm, double* vec){
    double sum;
    int i;
#pragma omp parallel private(sum) num_threads(4) if (ENABLE_PARALLEL)
    {
#pragma omp for nowait
        for (i = 0; i < crsm.nRows; i++){
            sum = 0.0;
            for (int j = crsm.rowindex[i]; j < crsm.rowindex[i + 1]; j++){
                sum += crsm.Value[j] * vec[crsm.col[j]];
            }
            (*result)[i] = sum;
        }
    }
}



   
