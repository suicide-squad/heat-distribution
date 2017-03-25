
#define _CRT_SECURE_NO_WARNINGS


#include <stdio.h>
#include <stdlib.h>
#include "crs.h"
#include <math.h>
#include <omp.h>

int nX;
void createCRS(CRSMatrix *crsm, double coeff1, double coeff2);

int main() {

    double xStart = 0.0, xEnd = 0.0;
    double sigma = 0.0;
    double tStart = 0.0, tFinal = 0.0;
    double dt = 0.0;
    int check = 0;
    volatile double *U = NULL;
    double stepp;
    int sizeTime = 0;



    FILE *fp;
    if ((fp = fopen("C:\\Users\\ролд\\Documents\\Visual Studio 2013\\Projects\\Eiler\\Eiler\\INPUT.txt", "r")) == NULL) {
        printf("File cant open!\n");
        exit(-1);
    }

    // Заполнение конфигурационных настроек
    fscanf(fp, "XSTART=%lf\n", &xStart);
    fscanf(fp, "XEND=%lf\n", &xEnd);
    fscanf(fp, "SIGMA=%lf\n", &sigma);
    fscanf(fp, "NX=%d\n", &nX);
    fscanf(fp, "TSTART=%lf\n", &tStart);
    fscanf(fp, "TFINISH=%lf\n", &tFinal);
    fscanf(fp, "dt=%lf\n", &dt);
    fscanf(fp, "BC=%d\n", &check);

    U = (double*)malloc(sizeof(double) * (nX + 2));

    // Заполнение функции в нулевой момент времени
    for (int i = 1; i < nX + 1; i++) {
        fscanf(fp, "%lf", &U[i]);
    }
    U[0] = U[nX + 1] = 0.0;

    fclose(fp);
    printf("%lf; %lf; %lf; %d; %lf; %lf; %lf; %d;\n", xStart, xEnd, sigma, nX, tStart, tFinal, dt, check);
    //  Вычисление шага по x
    stepp = fabs(xStart - xEnd) / (double)nX;

    if (2 * dt*sigma > stepp*stepp) {
        printf("Not correct stepp!");
        exit(-1);
    }
    sizeTime = (int)((tFinal - tStart) / dt);
    double backstepp = 1.0 / (stepp*stepp);
    CRSMatrix A;
    double coeff1 = backstepp;
    double coeff2 = -2.0*coeff1;
    createCRS(&A, coeff1, coeff2);

    CRSMatrix B;
    coeff1 = dt*backstepp*0.5;
    coeff2 = 1.0-2.0*coeff1;
    createCRS(&B, coeff1, coeff2);
    
    CRSMatrix C;
    coeff1 = dt*backstepp;
    coeff2 = 1.0 -2.0*coeff1;
    createCRS(&C, coeff1, coeff2);

    double* k1 = (double *)malloc((nX + 2)*sizeof(double));
    double* k2 = (double *)malloc((nX + 2)*sizeof(double));
    double* k3 = (double *)malloc((nX + 2)*sizeof(double));
    double* k4 = (double *)malloc((nX + 2)*sizeof(double));
    double* Unext = (double*)malloc(sizeof(double) * (nX + 2));
    double* tmp;
    double h = dt / 6.0;
    printf("start\n");

    double t0 = omp_get_wtime();
    for (int i = 1; i <= sizeTime; i++) {
        multCRSMatrix(&k1, A, U);
        multCRSMatrix(&k2, B, k1);
        multCRSMatrix(&k3, B, k2);
        multCRSMatrix(&k4, C, k3);
        int g;

#pragma omp parallel for num_threads(4) if (ENABLE_PARALLEL)
        for (g = 0; g < nX+2; g++)
            Unext[g] = U[g] + h*(k1[g] + 2.0*k2[g] + 2.0*k3[g] + k4[g]);        
            
        
        tmp = U;
        U = Unext;
        Unext = tmp;
    }
    double t1 = omp_get_wtime();
    double diffTime = t1 - t0;
    printf("Time\t%.15lf\n", diffTime);
    
    fp = 0;
    fp = fopen("C:\\Users\\ролд\\Desktop\\heat-distribution\\result\\svetaRungeSparse.txt", "w");


    if (!fp)
        return -1;


    for (int x = 1; x < nX + 1; x++)
        fprintf(fp, "%.15le\n", U[x]);
    printf("vse");

    fclose(fp);

    free(Unext);
    free(U);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    return 0;
}


void createCRS(CRSMatrix *crsm, double coeff1, double coeff2){
    initCRSMartix(nX + 2, 3 * nX + 2, crsm);
    crsm->Value[0] = coeff1;
    crsm->Value[1] = coeff2;
    crsm->Value[2] = coeff1;
    for (int i = 3; i < 3 * nX + 1; i += 3){
        crsm->Value[i] = coeff1;
        crsm->Value[i + 1] = coeff2;
        crsm->Value[i + 2] = coeff1;
    }
    crsm->Value[3 * nX + 1] = 1;
    crsm->col[0] = 0;
    crsm->col[1] = 1;
    crsm->col[2] = 2;
    crsm->col[3 * nX + 1] = nX + 1;
    int j = 0;
    for (int i = 3; i < 3 * nX + 1; i += 3){
        crsm->col[i] = j;
        crsm->col[i + 1] = j + 1;
        crsm->col[i + 2] = j + 2;
        j++;
    }
    crsm->rowindex[0] = 0;
    crsm->rowindex[1] = 1;
    for (int i = 2; i < nX + 2; i++)
        crsm->rowindex[i] = crsm->rowindex[i - 1] + 3;
    crsm->rowindex[nX + 2] = crsm->rowindex[nX + 1] + 1;

}