
#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"
#include "stdlib.h"
#include "math.h"




int main() {
    FILE *fp;

    int i, j;
    int sizeTime = 0;
    int curTime, prevTime;

    const int N = 2;
    double U[N];

    double xStart = 0.0, xEnd = 0.0;
    double sigma = 0.0;
    int nX = 1;
    double tStart = 0.0, tFinal = 0.0;
    double dt = 0.0;
    unsigned char check = 0;
    double step;


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
    fscanf(fp, "BC=%hhu\n", &check);

    nX = 500;
    sizeTime = (int)((tFinal - tStart) / dt);

    for (int x = 0; x < N; x++) {
        U[x] = (double *)malloc((nX + 2)*sizeof(double));
    }
    //  Вычисление шага по x
    step = fabs(xStart - xEnd) / nX;

    if (2 * dt*sigma > step*step) {
        printf("Not correct step!");
        exit(-1);
    }
    double backstep = 1 / (step*step);
    double h = dt/6;
    double polh = dt / 2;
    // Заполнение функции в нулевой момент времени
    for (int x = 1; x < nX - 1; x++) {
        fscanf(fp, "%lf", &U[0][x]);
    }

    fclose(fp);


    // Задание граничных условий
    if (check == 2) {
        U[0][0] = U[0][1];
        U[0][nX - 1] = U[0][nX - 2];
    }
    else if (check == 1) {
        printf("wtf");
    }

    for (int t = 1; t <= sizeTime; t++) {
        curTime = t%N;
        prevTime = (t + 1) % N;


    double* k1= (double *)malloc((nX + 2)*sizeof(double));
    for (int x = 1; x < nX + 2; x++){
        k1[x] = (U[prevTime][x - 1] + 2 * U[prevTime][x] + U[prevTime][x + 1])*backstep;
    }

    double* k2 = (double *)malloc((nX + 2)*sizeof(double));
    for (int x = 1; x < nX + 2; x++){
        k2[x] = (U[prevTime][x - 1] + k1[x - 1] * polh + 2 * U[prevTime][x] + k1[x] * dt + U[prevTime][x + 1] + k1[x + 1] * polh)*backstep;
    }

    double* k3 = (double *)malloc((nX + 2)*sizeof(double));
    for (int x = 1; x < nX + 2; x++){
        k3[x] = (U[prevTime][x - 1] + k2[x - 1] * polh + 2 * U[prevTime][x] + k2[x] * dt + U[prevTime][x + 1] + k2[x + 1] * polh)*backstep;
    }

    double* k4 = (double *)malloc((nX + 2)*sizeof(double));
    for (int x = 1; x < nX + 2; x++){
        k4[x] = (U[prevTime][x - 1] + k3[x - 1] * polh + 2 * U[prevTime][x] + k3[x] * dt + U[prevTime][x + 1] + k3[x + 1] * polh)*backstep;
    }
   

    for (int x = 1; x < nX + 2; x++){
        U[curTime][x] = U[prevTime][x] + h*(k1[x] + 2 * k2[x] + 2 * k3[x] + k4[x]);
    }
        //граничные условия
        if (check == 2) {
            U[curTime][0] = U[curTime][1];
            U[curTime][nX - 1] = U[curTime][nX - 2];
        }
        else if (check == 1) {
            printf("wtf");
        }
    }
        fp = fopen("C:\\Users\\ролд\\Documents\\Visual Studio 2013\\Projects\Eiler\\SvetaRungeKutti\\svetaRK.txt", "w");
        for (int x = 1; x < nX + 1; x++)
            fprintf(fp, "%.15le\n", U[sizeTime % 2][x]);
        fclose(fp);
        for (int i = 0; i < 2; i++)
            free(U[i]);

        return 0;
    }