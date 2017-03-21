//
// Created by kirill on 24.10.16.
//


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <sp_mat.h>

const double EPS = 1e-10;

const char pathInput[]  = "../../../../../../initial/INPUT.txt";
const char pathResult[] = "../../../../../../result/Kirill/ImplicitEuler.txt";

int init(double *, double *, double *, double *, double *, double *, int *, TYPE **);
void createSpMat(spMatrix *, TYPE);
bool dist(double *, double *, size_t N);

void final(TYPE *, const char *path);

size_t nX;

int main(int argc, char** argv) {
  double xStart, xEnd;
  double sigma;
  double tStart, tFinal;
  double dt;
  int check;

  TYPE* U;

  //------------------------------------------------------------------------
  //                      Инициализация данных
  //------------------------------------------------------------------------

  if ( init(&xStart, &xEnd, &sigma, &tStart, &tFinal, &dt, &check, &U) )
    return -1;

  double step = fabs(xStart - xEnd) / nX;
  size_t sizeTime = (size_t)((tFinal - tStart) / dt);

  #if ENABLE_PARALLEL
    printf("ПАРАЛЛЕЛЬНАЯ ВЕРСИЯ!\n");
  #endif
  printf("TIMESIZE = %lu; NX = %lu\n", sizeTime, nX);

  //------------------------------------------------------------------------
  //              Заполнение значений и номера столбцов матрицы
  //------------------------------------------------------------------------

  spMatrix A;

  double coeff = sigma*dt/(step*step);
  createSpMat(&A, coeff);

  // -----------------------------------------------------------------------
  //                              Вычисления
  //------------------------------------------------------------------------

  TYPE* X1 = (TYPE*)malloc(sizeof(TYPE) * (nX + 2));
  TYPE* X2 = (TYPE*)malloc(sizeof(TYPE) * (nX + 2));

  TYPE* tmp;

  double k = (step*step)/(step*step + 2*sigma*dt);

  int count = 0;
  double t0 = omp_get_wtime();

  for (int i = 1; i <= sizeTime; i++) {
    memcpy(X1, U, (nX + 2)*sizeof(double));

    do {
      multMV(&X2, A, X1);

      for (int j = 1; j < nX + 1; j++) {
        X2[j] = (U[j] + X2[j])*k;
      }

      X2[0] = X2[1];
      X2[nX + 1] = X2[nX];

      tmp = X1;
      X1 = X2;
      X2 = tmp;

      count++;

    } while (!dist(X1, X2, nX + 2));

    memcpy(U, X1, (nX + 2)*sizeof(double));

  }


  double t1 = omp_get_wtime();
  printf("\nfinish!\n\n");

  //------------------------------------------------------------------------
  //                       Вывод результатов и чистка памяти
  //------------------------------------------------------------------------

  double diffTime = t1 - t0;
  double gflop = (2*2*nX)*count*sizeTime*1.0/1000000000.0;

  printf("EPS\t%.le\n", EPS);
  printf("COUNT\t%d\n", count);
  printf("TIME\t%.15lf\n", diffTime);
//  printf("GFlop\t%.lf\n", gflop);
//  printf("GFlop's\t%.15lf\n", gflop*1.0/diffTime);


  final(U, pathResult);

  free(U);
  free(X1);
  free(X2);

  freeSpMat(&A);
  return 0;

}

/*
____________________________________________________________________________


                          РЕАЛИЗАЦИЯ ФУНКЦИЙ

 ____________________________________________________________________________

*/


int init(double *xStart, double *xEnd, double *sigma, double *tStart, double *tFinal, double *dt, int *check, TYPE **U) {
  FILE *fp;
  if ((fp = fopen(pathInput, "r")) == NULL) {
    printf("Не могу найти файл!\n");
    return -2;
  };

  if ( !fscanf(fp, "XSTART=%lf\n", xStart) )
    return -1;
  if ( !fscanf(fp, "XEND=%lf\n", xEnd) )
    return -1;
  if ( !fscanf(fp, "SIGMA=%lf\n", sigma) )
    return -1;
  if ( !fscanf(fp, "NX=%lu\n", &nX) )
    return -1;
  if ( !fscanf(fp, "TSTART=%lf\n", tStart) )
    return -1;
  if ( !fscanf(fp, "TFINISH=%lf\n", tFinal) )
    return -1;
  if ( !fscanf(fp, "dt=%lf\n", dt) )
    return -1;
  if ( !fscanf(fp, "BC=%d\n", check) )
    return -1;

  *U = (TYPE*)malloc(sizeof(TYPE) * (nX + 2));

  // Заполнение функции в нулевой момент времени
  for (int i = 1; i < nX - 1; i++)
    if ( !fscanf(fp, "%lf", &(*U)[i]) )
      return -1;
  (*U)[0] = (*U)[nX + 1] = 0.0;
  fclose(fp);

  return 0;
}

void createSpMat(spMatrix *mat, TYPE coeff) {

  size_t nz = nX*2;
  size_t size = nX + 2;
  initSpMat(mat, nz, size);

  int j = 0;
  for (int i = 0; i < nz; i += 2) {
    mat->value[i] = coeff;       mat->col[i] = j;
    mat->value[i + 1] = coeff;   mat->col[i + 1] = j + 2;
    j++;
  }

  mat->rowIndex[0] = 0;
  mat->rowIndex[1] = 3;
  for (int i = 2; i < size; i++) {
    mat->rowIndex[i] = mat->rowIndex[i - 1] + 2;
  }
  mat->rowIndex[size] = mat->rowIndex[size - 1];
}

void final(TYPE *UFin, const char *path) {
  FILE *fp;
  fp = fopen(path, "w");

  for (int i = 1; i < nX + 1; i++)
    fprintf(fp, "%.15le\n", UFin[i]);

  fclose(fp);
}

bool dist(double *x1, double *x2, size_t N) {
  for (int i = 0; i < N; i++)
    if (fabs(x1[i] - x2[i]) > EPS)
      return false;
  return true;
}

