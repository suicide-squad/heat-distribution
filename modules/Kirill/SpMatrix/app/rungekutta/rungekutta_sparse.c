//
// Created by kirill on 24.10.16.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "sp_mat.h"

#define ROOT 0

const char pathInput[]  = "../../../../../../initial/INPUT.txt";
const char pathResult[] = "../../../../../../result/Kirill/runge1D_MPI.txt";

int init(double *, double *, double *, size_t *nX, double *, double *, double *, int *, TYPE **);
void createSpMat(spMatrix *, int nX, int reserve, TYPE, TYPE);
void final(TYPE *, size_t nX, const char *path);

size_t nX;

int main(int argc, char **argv) {
  double t0 = 0.0, t1 = 0.0;
  size_t nX, sizeTime;

  double step, dt;
  TYPE* U = NULL, *Ur = NULL;

  MPI_Status status;
  int sizeP, rankP;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankP);
  //------------------------------------------------------------------------
  //                       Инициализация данных
  //------------------------------------------------------------------------

  if (rankP == ROOT) {
    double xStart, xEnd;
    double sigma;
    double tStart, tFinal;
    int check;
    if (init(&xStart, &xEnd, &sigma, &nX, &tStart, &tFinal, &dt, &check, &U))
      return -1;

    step = fabs(xStart - xEnd) / nX;
    sizeTime = (size_t) ((tFinal - tStart) / dt);
    if (2 * sigma * dt > step * step) {
      printf("Выбор шага по времени не возможен в силу условия устойчивости!\n");
      printf("%.10lf > %.10lf\n", 2 * sigma * dt, step * step);
      printf("Предлагаю взять dt = %.10lf\n", step * step / (2.0 * sigma));
      return -1;
    }
    #if ENABLE_PARALLEL
      printf("ПАРАЛЛЕЛЬНАЯ ВЕРСИЯ!\n");
    #endif
    printf("TIMESIZE = %lu; NX = %lu\n", sizeTime, nX);
  }

  MPI_Bcast(&sizeTime, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&nX, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&step, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&dt, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

  int block = (int)nX/sizeP;

  int reserve = 2;
  int shift = 1;
  if (sizeP != 1) {
    if (rankP == ROOT) {
      // краевое условие + 4 элемента
      reserve = 5;
      shift = 1;
    }
    else if ( rankP == sizeP - 1 ){
    // 4 элемента + 4 элемента
      reserve = 5;
      shift = 4;
    } else {
      reserve = 8;
      shift = 4;
    }
  }

  Ur = (TYPE*)malloc(sizeof(TYPE)*(block + reserve));
  MPI_Scatter(U, block, MPI_DOUBLE, Ur + shift, block, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  Ur[0] = Ur[1];
  Ur[block + reserve - 1] = Ur[block];

  //------------------------------------------------------------------------
  //          Заполнение значений и номера столбцов матрицы
  //------------------------------------------------------------------------

  spMatrix A;
  double coeff1 = 1.0/(step*step);
  double coeff2 = -2.0*coeff1;
  initSpMat(&A, (block + reserve)*3, block + reserve);
  createSpMat(&A, block, reserve, coeff1, coeff2);

  spMatrix B;
  coeff1 = dt*coeff1*0.5;
  coeff2 = 1.0 - 2.0*coeff1;
  initSpMat(&B, (block + reserve)*3, block + reserve);
  createSpMat(&B, block, reserve, coeff1, coeff2);

  spMatrix C;
  coeff1 = coeff1*2.0;
  coeff2 = 1.0 - 2.0*coeff1;
  initSpMat(&C, (block + reserve)*3, block + reserve);
  createSpMat(&C, block, reserve, coeff1, coeff2);

  // -----------------------------------------------------------------------
  //                         Вычисления
  //------------------------------------------------------------------------

  TYPE* UrNext = (TYPE*)malloc(sizeof(TYPE) * (block + reserve));
  TYPE* k1 = (TYPE*)malloc(sizeof(TYPE) * (block + reserve));
  TYPE* k2 = (TYPE*)malloc(sizeof(TYPE) * (block + reserve));
  TYPE* k3 = (TYPE*)malloc(sizeof(TYPE) * (block + reserve));
  TYPE* k4 = (TYPE*)malloc(sizeof(TYPE) * (block + reserve));

  TYPE* tmp;

  double h = dt/6.0;

  if (rankP == ROOT) t0 = omp_get_wtime();
  for (int i = 1; i <= sizeTime; i++) {

    if (sizeP != 1) {
      if (rankP == ROOT) {
        MPI_Send(Ur + block + reserve - 8, 4, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(Ur + block + reserve - 4, 4, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD, &status);
      } else if (rankP == sizeP - 1) {
        MPI_Send(Ur + 4, 4, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(Ur + 0, 4, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, &status);
      } else {
        MPI_Send(Ur + block + reserve - 8, 4, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(Ur + block + reserve - 4, 4, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(Ur + 4, 4, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(Ur + 0, 4, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, &status);
      }
    }

    // k1 = A*U
    multMV(&k1, A, Ur);

    // k2 = B*k1
    multMV(&k2, B, k1);

    // k3 = B*k2
    multMV(&k3, B, k2);

    // k4 = C*k3
    multMV(&k4, C, k3);

    // UrNext = U + (k1 + k2*2 + k3*2 + k4)*h;
    sumV(block + reserve, h, &UrNext, Ur, k1, k2, k3, k4);

    tmp = Ur;
    Ur = UrNext;
    UrNext = tmp;
  }
  if (rankP == ROOT) {
    t1 = omp_get_wtime();
    printf("\nfinish!\n\n");
  }

  //------------------------------------------------------------------------
  //                       Вывод результатов и чистка памяти
  //------------------------------------------------------------------------

  //  Сбор информацию в итоговый массив
  MPI_Gather(Ur + shift, block, MPI_DOUBLE, U, block, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

  if (rankP == ROOT) {
    double diffTime = t1 - t0;
    double gflop = (2 * 3 * nX * 4 + 7 * (nX + 2)) * sizeTime * 1.0 / 1000000000.0;
    printf("Time\t%.15lf\n", diffTime);
    printf("GFlop\t%.lf\n", gflop);
    printf("GFlop's\t%.15lf\n", gflop * 1.0 / diffTime);

    final(U, nX, pathResult);
    free(U);
  }
  free(UrNext);
  free(k1);
  free(k2);
  free(k3);
  free(k4);

  freeSpMat(&A);
  freeSpMat(&B);
  freeSpMat(&C);

  MPI_Finalize();
  return 0;

}

/*
____________________________________________________________________________

                          РЕАЛИЗАЦИЯ ФУНКЦИЙ
 ____________________________________________________________________________

*/


int init(double *xStart,
         double *xEnd,
         double *sigma,
         size_t *nX,
         double *tStart,
         double *tFinal,
         double *dt,
         int *check,
         TYPE **U) {
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
  if ( !fscanf(fp, "NX=%lu\n", nX) )
    return -1;
  if ( !fscanf(fp, "TSTART=%lf\n", tStart) )
    return -1;
  if ( !fscanf(fp, "TFINISH=%lf\n", tFinal) )
    return -1;
  if ( !fscanf(fp, "dt=%lf\n", dt) )
    return -1;
  if ( !fscanf(fp, "BC=%d\n", check) )
    return -1;

  *U = (TYPE*)malloc(sizeof(TYPE) * (*nX));

  // Заполнение функции в нулевой момент времени
  for (int i = 0; i < *nX; i++)
    if ( !fscanf(fp, "%lf", &(*U)[i]) )
      return -1;
  fclose(fp);

  return 0;
}

void createSpMat(spMatrix *mat, int nX, int reserve, TYPE coeff1, TYPE coeff2) {
  int j = 0;
  mat->value[0] = coeff1;  mat->col[0] = 0;
  mat->value[1] = coeff2; mat->col[1] = 1;
  mat->value[2] = coeff1;  mat->col[2] = 2;
  for (int i = 3; i < (nX+reserve)*3 - 3; i += 3) {
    mat->value[i] = coeff1;       mat->col[i] = j++;
    mat->value[i + 1] = coeff2;  mat->col[i + 1] = j++;
    mat->value[i + 2] = coeff1;   mat->col[i + 2] = j--;
  }
  mat->value[(nX+reserve)*3 - 3] = coeff1; mat->col[(nX+reserve)*3 - 3] = nX + reserve - 3;
  mat->value[(nX+reserve)*3 - 2] = coeff2; mat->col[(nX+reserve)*3 - 2] = nX + reserve - 2;
  mat->value[(nX+reserve)*3 - 1] = coeff1; mat->col[(nX+reserve)*3 - 1] = nX + reserve - 1;

  mat->rowIndex[0] = 0;
  for (int i = 1; i < nX + reserve; i++) {
    mat->rowIndex[i] = mat->rowIndex[i - 1] + 3;
  }
  mat->rowIndex[nX + reserve] = mat->rowIndex[nX + reserve - 1] + 3;
}

void final(TYPE *UFin, size_t nX, const char *path) {
  FILE *fp;
  fp = fopen(path, "w");

  for (int i = 0; i < nX; i++)
    fprintf(fp, "%.15le\n", UFin[i]);

  fclose(fp);
}
