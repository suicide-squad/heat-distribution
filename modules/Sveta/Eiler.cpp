
#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <omp.h>

const int ENABLE_PARALLEL = 1;

int main() {
    FILE *fp;

    int i, j;
    int sizeTime = 0;
    int curTime, prevTime;

    const int N = 2;
    double* U[N];

    double xStart = 0.0, xEnd = 0.0;
    double sigma = 0.0;
    int nX = 1;
    double tStart = 0.0, tFinal = 0.0;
    double dt = 0.0;
    int check = 0;
    double step = 0.0;

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

    sizeTime = (int)((tFinal - tStart) / dt);

    
    for (i = 0; i < N; i++) {
        U[i] = (double *)malloc((nX + 2)*sizeof(double));
    }

    //  Вычисление шага по x
    step = fabs(xStart - xEnd) / nX;

    if (2 * dt*sigma > step*step) {
        printf("Not correct step!");
        exit(-1);
    }
    int backstep = 1 / (step*step);

    // Заполнение функции в нулевой момент времени
    for (i = 1; i < nX + 1; i++) {
        fscanf(fp, "%lf", &U[0][i]);
    }

    fclose(fp);

    // Задание граничных условий
    if (check == 2) {
        U[0][0] = U[0][1];
        U[0][nX + 1] = U[0][nX];
    }
    else if (check == 1) {
        // ???
        printf("HZ");
    }
    printf("%d, %d;\n", sizeTime, nX);
    double coef = dt*sigma*backstep;
    // Заполнение сетки
    double t0 = omp_get_wtime();
    for (i = 1; i <= sizeTime; i++) {
        curTime = i%N;
        prevTime = (i + 1) % N;
        omp_set_num_threads(4);
        {

#pragma omp parallel for if (ENABLE_PARALLEL)
            for (j = 1; j < nX + 1; j++)
                U[curTime][j] = coef*(U[prevTime][j - 1] - 2 * U[prevTime][j] + U[prevTime][j + 1]) + U[prevTime][j];
        }
        // Задание граничных условий
        if (check == 2) {
            U[curTime][0] = U[curTime][1];
            U[curTime][nX+1] = U[curTime][nX ];
        }
        else if (check == 1) {
            // ???
            
        }

    }
    double t1 = omp_get_wtime();
    double diffTime = t1 - t0;


    printf("Time\t%.15lf\n", diffTime);
    fp = fopen("C:\\Users\\ролд\\Documents\\Visual Studio 2013\\Projects\\Eiler\\Eiler\\svetaEuler.txt", "w");

    // Вывод результатов
    for (j = 1; j < nX + 1; j++)
        fprintf(fp, "%.15le\n", U[sizeTime % 2][j]);

    fclose(fp);

    // Чистка мусора
    
    
    free(U);
    return 0;
}