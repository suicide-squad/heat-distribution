//
// Created by kirill on 24.10.16.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sp_mat.h>

int init(double *, double *, double *, double *, double *, double *, int *, TYPE **);
void createSpMat(spMatrix *, TYPE);
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
  //                      Инициализация данных
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
  //              Заполнение значений и номера столбцов матрицы
  //------------------------------------------------------------------------

  spMatrix A;
  double coeff1 = dt/(step*step);
  double coeff2 = 1 - 2.0*coeff1;
  createSpMat(&A, coeff1);

  // -----------------------------------------------------------------------
  //                              Вычисления
  //------------------------------------------------------------------------

  TYPE* UNext = (TYPE*)malloc(sizeof(TYPE) * (nX + 2));
  TYPE* tmp;

  double t0 = omp_get_wtime();
  for (int i = 1; i <= sizeTime; i++) {
    // UNext = A*U
    multMV(&UNext, A, U);

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
  unsigned long long flop = (2*3*nX + 2*2)*sizeTime;
  printf("Time\t%.15lf\n", diffTime);
  printf("Flop\t%.0llu\n", flop);
  printf("GFlops\t%.15lf\n", flop*1.0/(diffTime*1000000000.0));

  final(&U);

//  free(U);
  free(UNext);

//  freeSpMat(&A);
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

void createSpMat(spMatrix *mat, TYPE coeff) {

  initSpMat(mat, nX*3, nX + 3);

  int j = 0;
    mat->value[0] = 1.0;          mat->col[0] = 0;
  for (int i = 1; i < 3*nX + 1; i += 3) {
    mat->value[i] = coeff;       mat->col[i] = j++;
    mat->value[i + 1] = coeff2;   mat->col[i + 1] = j++;
    mat->value[i + 2] = coeff;   mat->col[i + 2] = j--;
  }
    mat->value[nX*3 + 1] = 1.0;   mat->col[3*nX + 1]  = (int)nX + 1;

  mat->rowIndex[0] = 0;
  mat->rowIndex[1] = 1;
  for (int i = 2; i < nX + 2; i++) {
    mat->rowIndex[i] = mat->rowIndex[i - 1] + 3;
  }
  mat->rowIndex[nX + 2] = mat->rowIndex[nX + 1] + 1;
}

int final(TYPE *UFin) {
  FILE *fp;
  fp = fopen("./../../../../result/kirillEulerSparse.txt", "w");

  for (int i = 1; i < nX + 1; i++)
    fprintf(fp, "%.15le\n", (*UFin)[i]);

  fclose(fp);
}
