//
// Created by kirill on 14.11.16.
//

#include <stdio.h>
#include <math.h>

#include <stdlib.h>
#include <sp_mat.h>

#ifdef _COMPLEX_

int init(double *, double *, double *, double *, double *, double *, int *, TYPE **);
void createSpMat(spMatrix *, TYPE);

void printInFile(TYPE *UFin, FILE* fpRe, FILE* fpIm);
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

  printf("TIMESIZE = %lu; NX = %lu\n", sizeTime, nX);

  //------------------------------------------------------------------------
  //              Заполнение значений и номера столбцов матрицы
  //------------------------------------------------------------------------

  spMatrix A1;
  complex k1exp = CMPLX(0., 1.);

  printf("%lf\t%+.lfi\n", creal(k1exp), cimag(k1exp));

  createSpMat(&A1, k1exp);

  spMatrix A2;

  complex k2exp = CMPLX(1., dt*0.5);

  printf("%.15lf\t%+.15lfi\n", creal(k2exp), cimag(k2exp));

  createSpMat(&A2, k2exp);

  spMatrix A3;

//  complex k3exp = CMPLX(0., dt*0.5) + 1./k2exp;
complex k3exp = k2exp;

  printf("%.15lf\t%+.15lfi\n", creal(k3exp), cimag(k3exp));

  createSpMat(&A3, k3exp);

  spMatrix A4;

//  complex k4exp = CMPLX(0., dt) + 1./(k2exp*k3exp);
complex k4exp = CMPLX(1., dt);
  printf("%.15lf\t%+.15lfi\n", creal(k4exp), cimag(k4exp));

  createSpMat(&A4, k4exp);

  FILE* fpRe = fopen("./../../../../result/complex/Re.txt", "w");
  FILE* fpIm = fopen("./../../../../result/complex/Im.txt", "w");

  // -----------------------------------------------------------------------
  //                              Вычисления
  //------------------------------------------------------------------------

  TYPE* UNext = (TYPE*)malloc(sizeof(TYPE) * nX);

  TYPE* k1 = (TYPE*)malloc(sizeof(TYPE) * nX);
  TYPE* k2 = (TYPE*)malloc(sizeof(TYPE) * nX);
  TYPE* k3 = (TYPE*)malloc(sizeof(TYPE) * nX);
  TYPE* k4 = (TYPE*)malloc(sizeof(TYPE) * nX);

  TYPE* tmp;

  printInFile(U, fpRe, fpIm);

  double t0 = omp_get_wtime();
  for (int i = 1; i <= sizeTime; i++) {
    // k1 = A1*U
    multMV(&k1, A1, U);

    // k2 = A2*k1
    multMV(&k2, A2, k1);

    // k3 = A3*k2
    multMV(&k3, A3, k2);

    // k4 = A4*k3
    multMV(&k4, A4, k3);

    tmp = U;
    U = UNext;
    UNext = tmp;

    if ( i%10000 == 0 )
      printInFile(U, fpRe, fpIm);
  }
  double t1 = omp_get_wtime();

  fclose(fpRe);
  fclose(fpIm);

  printf("finish!\n\n");

  //------------------------------------------------------------------------
  //                       Вывод результатов и чистка памяти
  //------------------------------------------------------------------------

  double diffTime = t1 - t0;
  printf("Time\t%.15lf\n", diffTime);

  final(U);

  free(U);
  free(UNext);

  freeSpMat(&A1);
  freeSpMat(&A2);
  freeSpMat(&A3);
  freeSpMat(&A4);
  return 0;

}

/*
____________________________________________________________________________

                          РЕАЛИЗАЦИЯ ФУНКЦИЙ
 ____________________________________________________________________________

*/


int init(double *xStart, double *xEnd, double *sigma, double *tStart,
         double *tFinal, double *dt, int *check, TYPE **U) {
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

  *U = (TYPE *) malloc(sizeof(TYPE) * nX);

  double re, im = 0.;

  // Заполнение функции в нулевой момент времени
  for(int i = 0; i < nX; i++) {
    fscanf(fp, "%lf", &re);
    (*U)[i] = CMPLX(re, 0);
  }
  fclose(fp);

  return 0;
}
void createSpMat(spMatrix *mat, TYPE coeff) {

  initSpMat(mat, nX, nX);

  int j = 0;
  for (int i = 0; i < nX; i ++) {
    mat->value[i] = coeff;
    mat->col[i] = j++;
  }

  mat->rowIndex[0] = 0;
  for (int i = 1; i < nX + 1; i++) {
    mat->rowIndex[i] = i;
  }

}

void printInFile(TYPE *UFin, FILE* fpRe, FILE* fpIm) {
  for (int i = 0; i < nX; i++) {
    fprintf(fpRe, "%lf ", creal(UFin[i]));
    fprintf(fpIm, "%lf ", cimag(UFin[i]));
  }
  fprintf(fpRe, "\n");
  fprintf(fpIm, "\n ");
}

int final(TYPE *UFin) {
  FILE *fp;
  fp = fopen("./../../../../result/complex/RungeSparse.txt", "w");

  for (int i = 0; i < nX; i++)
    fprintf(fp, "%.15le\t%+.15lei\n", creal(UFin[i]), cimag(UFin[i]));

  fclose(fp);
}

#endif