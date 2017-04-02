//
// Created by kirill on 24.10.16.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sp_mat.h>
#include <mpi.h>
#include <string.h>

#define ROOT 0

const char pathInput[]  = "../../../../../../initial/INPUT.txt";
const char pathResult[] = "../../../../../../result/Kirill/euler1D_MPI.txt";

int init(double *, double *, double *, size_t *nX, double *, double *, double *, int *, TYPE **);
void createSpMat(spMatrix *, int nX, TYPE, TYPE);
void final(TYPE *, size_t nX, const char *path);


int main(int argc, char **argv) {
  size_t nX, sizeTime;

  double step, dt;
  TYPE* U = NULL, *Ur = NULL;

  MPI_Status status;
  int sizeP, rankP;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

  //------------------------------------------------------------------------
  //                      Инициализация данных
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

  Ur = (TYPE*)malloc(sizeof(TYPE)*(block + 2));
  MPI_Scatter(U, block, MPI_DOUBLE, Ur + 1, block, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  Ur[0] = Ur[1];
  Ur[block + 1] = Ur[block];

  //------------------------------------------------------------------------
  //              Заполнение значений и номера столбцов матрицы
  //------------------------------------------------------------------------

  spMatrix A;
  initSpMat(&A, (block + 2)*3, block + 2);

  double coeff1 = dt/(step*step);
  double coeff2 = 1.0 - 2.0*coeff1;
  createSpMat(&A, block, coeff1, coeff2);

  // -----------------------------------------------------------------------
  //                              Вычисления
  //------------------------------------------------------------------------

  TYPE* UrNext = (TYPE*)malloc(sizeof(TYPE) * (block + 2));
  TYPE* tmp;

  double t0 = omp_get_wtime();
  for (int i = 1; i <= sizeTime; i++) {
    // Передача соседям информацию о границах
    if (sizeP != 1) {
      if (rankP == ROOT) {
        MPI_Send(Ur + block + 0, 1, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(Ur + block + 1, 1, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD, &status);
      } else if (rankP == sizeP - 1) {
        MPI_Send(Ur + 1, 1, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(Ur + 0, 1, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, &status);
      } else {
        MPI_Send(Ur + block + 0, 1, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(Ur + block + 1, 1, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(Ur + 1, 1, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(Ur + 0, 1, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, &status);
      }
    }

    // UNext(rank) = A(rank)*U(rank)
    multMV(&UrNext, A, Ur);

    tmp = Ur;
    Ur = UrNext;
    UrNext = tmp;
  }

  //  Сбор информацию в итоговый массив
  MPI_Gather(Ur + 1, block, MPI_DOUBLE, U, block, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

  if (rankP == ROOT) {
    double t1 = omp_get_wtime();
    printf("\nfinish!\n\n");

    //------------------------------------------------------------------------
    //                       Вывод результатов и чистка памяти
    //------------------------------------------------------------------------

    double diffTime = t1 - t0;
    double gflop = (2 * 3 * nX + 2 * 2) * sizeTime * 1.0 / 1000000000.0;
    printf("Time\t%.15lf\n", diffTime);
    printf("GFlop\t%.lf\n", gflop);
    printf("GFlop's\t%.15lf\n", gflop * 1.0 / diffTime);

    final(U, nX, pathResult);
    free(U);
  }

  free(Ur);
  free(UrNext);
  freeSpMat(&A);

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

void createSpMat(spMatrix *mat, int nX, TYPE coeff, TYPE coeff2) {
  int j = 0;
  mat->value[0] = coeff;  mat->col[0] = 0;
  mat->value[1] = coeff2; mat->col[1] = 1;
  mat->value[2] = coeff;  mat->col[2] = 2;
  for (int i = 3; i < (nX + 2)*3 - 3; i += 3) {
    mat->value[i] = coeff;       mat->col[i] = j++;
    mat->value[i + 1] = coeff2;  mat->col[i + 1] = j++;
    mat->value[i + 2] = coeff;   mat->col[i + 2] = j--;
  }
  mat->value[(nX+2)*3 - 3] = coeff;  mat->col[(nX+2)*3 - 3] = (int)nX - 1;
  mat->value[(nX+2)*3 - 2] = coeff2; mat->col[(nX+2)*3 - 2] = (int)nX;
  mat->value[(nX+2)*3 - 1] = coeff;  mat->col[(nX+2)*3 - 1] = (int)nX + 1;

  mat->rowIndex[0] = 0;
  for (int i = 1; i < nX + 2; i++) {
    mat->rowIndex[i] = mat->rowIndex[i - 1] + 3;
  }
  mat->rowIndex[nX + 2] = mat->rowIndex[nX + 1] + 3;
}

void final(TYPE *UFin, size_t nX, const char *path) {
  FILE *fp;
  fp = fopen(path, "w");

  for (int i = 0; i < nX; i++)
    fprintf(fp, "%.15le\n", UFin[i]);

  fclose(fp);
}
