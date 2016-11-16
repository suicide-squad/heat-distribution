//
// Created by kirill on 24.10.16.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sp_mat.h"

int init(double *, double *, double *, double *, double *, double *, int *, TYPE **);
void createSpMat(spMatrix *, TYPE, TYPE);
int final(TYPE *);

size_t nX;

int main() {
  double xStart, xEnd;
  double sigma;
  double tStart, tFinal;
  double dt;
  int check;

  TYPE* U;


  //------------------------------------------------------------------------
  //                       Инициализация данных
  //------------------------------------------------------------------------

  init(&xStart, &xEnd, &sigma, &tStart, &tFinal, &dt, &check, &U);

  double step = fabs(xStart - xEnd) / nX;
  size_t sizeTime = (size_t)((tFinal - tStart) / dt);
  if (2 * sigma * dt > step * step) {
    printf("Выбор шага по времени не возможен в силу условия устойчивости!\n");
    printf("%.10lf > %.10lf\n", 2 * sigma * dt, step * step);
    printf("Предлагаю взять dt = %.10lf\n", step * step / (2.0 * sigma));

    return -1;
  }

  printf("TIMESIZE = %lu; NX = %lu\n", sizeTime, nX);

  //------------------------------------------------------------------------
  //          Заполнение значений и номера столбцов матрицы
  //------------------------------------------------------------------------

  spMatrix A;
  double coeff1 = 1.0/(step*step);
  double coeff2 = -2.0*coeff1;
  createSpMat(&A, coeff1, coeff2);

  spMatrix B;
  coeff1 = dt*coeff1*0.5;
  coeff2 = 1.0 - 2.0*coeff1;
  createSpMat(&B, coeff1, coeff2);

  spMatrix C;
  coeff1 = coeff1*2.0;
  coeff2 = 1.0 - 2.0*coeff1;
  createSpMat(&C, coeff1, coeff2);

  // -----------------------------------------------------------------------
  //                         Вычисления
  //------------------------------------------------------------------------

  TYPE* UNext = (TYPE*)malloc(sizeof(TYPE) * (nX + 2));
  TYPE* k1 = (TYPE*)malloc(sizeof(TYPE) * (nX + 2));
  TYPE* k2 = (TYPE*)malloc(sizeof(TYPE) * (nX + 2));
  TYPE* k3 = (TYPE*)malloc(sizeof(TYPE) * (nX + 2));
  TYPE* k4 = (TYPE*)malloc(sizeof(TYPE) * (nX + 2));

  TYPE* tmp;

  double h = dt/6.0;

  double t0 = omp_get_wtime();
  for (int i = 1; i <= sizeTime; i++) {

    // k1 = A*U
    multMV(&k1, A, U);

    // k2 = B*k1
    multMV(&k2, B, k1);

    // k3 = B*k2
    multMV(&k3, B, k2);

    // k4 = C*k3
    multMV(&k4, C, k3);

    // UNext = U + (k1 + k2*2 + k3*2 + k4)*h;
    sumV(nX + 2, h, &UNext, U, k1, k2, k3, k4);

    tmp = U;
    U = UNext;
    UNext = tmp;
  }
  double t1 = omp_get_wtime();
  printf("finish!\n\n");

  //------------------------------------------------------------------------
  //                       Вывод результатов и чистка памяти
  //------------------------------------------------------------------------

  double diffTime = t1 - t0;
  unsigned long long flop = (2*3*nX*4 + 7*(nX+2))*sizeTime;
  printf("Time\t%.15lf\n", diffTime);
  printf("Flop\t%.0llu\n", flop);
  printf("GFlops\t%.15lf\n", flop*1.0/(diffTime*1000000000.0));

  final(U);

  free(U);
  free(UNext);
  free(k1);
  free(k2);
  free(k3);
  free(k4);

  freeSpMat(&A);
  freeSpMat(&B);
  freeSpMat(&C);

  return 0;

}

/*
____________________________________________________________________________

                          РЕАЛИЗАЦИЯ ФУНКЦИЙ
 ____________________________________________________________________________

*/


int init(double *xStart, double *xEnd, double *sigma, double *tStart, double *tFinal, double *dt, int *check, TYPE **U) {
  FILE *fp;
  if ((fp = fopen("./../../../../initial/INPUT.txt", "r")) == NULL) {
    printf("Не могу найти файл!\n");
    return -1;
  };

  fscanf(fp, "XSTART=%lf\n", xStart);
  fscanf(fp, "XEND=%lf\n", xEnd);
  fscanf(fp, "SIGMA=%lf\n", sigma);
  fscanf(fp, "NX=%lu\n", &nX);
  fscanf(fp, "TSTART=%lf\n", tStart);
  fscanf(fp, "TFINISH=%lf\n", tFinal);
  fscanf(fp, "dt=%lf\n", dt);
  fscanf(fp, "BC=%d\n", check);

  *U = (TYPE*)malloc(sizeof(TYPE) * (nX + 2));

  // Заполнение функции в нулевой момент времени
  for (int i = 1; i < nX - 1; i++)
    fscanf(fp, "%lf", &(*U)[i]);
  (*U)[0] = (*U)[nX + 1] = 0.0;
  fclose(fp);

  return 0;
}

void createSpMat(spMatrix *mat, TYPE coeff1, TYPE coeff2) {

  initSpMat(mat, nX*3, nX + 3);

  int j = 0;
  for (int i = 0; i < 3*nX; i += 3) {
    mat->value[i] = coeff1;       mat->col[i] = j++;
    mat->value[i + 1] = coeff2;   mat->col[i + 1] = j++;
    mat->value[i + 2] = coeff1;   mat->col[i + 2] = j--;
  }

  mat->rowIndex[0] = 0;
  mat->rowIndex[1] = 0;
  for (int i = 2; i < nX + 2; i++) {
    mat->rowIndex[i] = mat->rowIndex[i - 1] + 3;
  }
  mat->rowIndex[nX + 2] = mat->rowIndex[nX + 1];
}

int final(TYPE *UFin) {
  FILE *fp;
  fp = fopen("./../../../../result/kirillRungeKuttaSparse.txt", "w");

  for (int i = 1; i < nX + 1; i++)
    fprintf(fp, "%.15le\n", UFin[i]);

  fclose(fp);
}
