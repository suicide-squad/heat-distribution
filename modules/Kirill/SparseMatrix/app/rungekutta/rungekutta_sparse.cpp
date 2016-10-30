//
// Created by kirill on 24.10.16.
//

#include <stdio.h>
#include <math.h>

#include "sparse_matrix.h"

int init(double&, double&, double&, size_t&, double&, double&, double&, int&, vector&);
int final(const size_t, vector&);

int main() {
  double xStart, xEnd;
  double sigma;
  size_t nX;
  double tStart, tFinal;
  double dt;
  int check;

  vector U;

  //------------------------------------------------------------------------
  //                       Инициализация данных
  //------------------------------------------------------------------------

  init(xStart, xEnd, sigma, nX, tStart, tFinal, dt, check, U);

  double step = fabs(xStart - xEnd) / nX;
  int sizeTime = static_cast<int>((tFinal - tStart) / dt);
  if (2 * sigma * dt > step * step) {
    printf("Выбор шага по времени не возможен в силу условия устойчивости!\n");
    printf("%.10lf > %.10lf\n", 2 * sigma * dt, step * step);
    printf("Предлагаю взять dt = %.10lf\n", step * step / (2.0 * sigma));

    return -1;
  }

  printf("TIMESIZE = %d; NX = %lu\n", sizeTime, nX);


  // Для матрицы А
  double* valA = new double[nX*3];
  int* colA = new int[nX*3];
  int* rowIA = new int[nX + 3];

  // Для матрицы B
  double* valB = new double[nX*3];
  int* colB = new int[nX*3];
  int* rowIB = new int[nX + 3];

  // Для матрицы C
  double* valC = new double[nX*3];
  int* colC = new int[nX*3];
  int* rowIC = new int[nX + 3];

  //------------------------------------------------------------------------
  //          Заполнение значений и номера столбцов матрицы
  //------------------------------------------------------------------------

  double h = 1.0/(step*step);
  double r1 = dt*h*0.5;
  double r2 = dt*h;

  //          Matrix A

  int j = 0;
  for (int i = 0; i < 3*nX; i +=3 ) {
    valA[i] = h;            colA[i] = j++;
    valA[i + 1] = -2.0*h;   colA[i + 1] = j++;
    valA[i + 2] = h;        colA[i + 2] = j--;
  }

  rowIA[0] = 0;
  rowIA[1] = 0;
  for (int i = 2; i < nX + 2; i++) {
    rowIA[i] = rowIA[i - 1] + 3;
  }
  rowIA[nX + 2] = rowIA[nX + 1];

  //          Matrix B
  j = 0;
  for (int i = 0; i < 3*nX; i += 3) {
    valB[i] = r1;                 colB[i] = j++;
    valB[i + 1] = 1.0 - 2.0*r1;   colB[i + 1] = j++;
    valB[i + 2] = r1;             colB[i + 2] = j--;
  }

  rowIB[0] = 0;
  rowIB[1] = 0;
  for (int i = 2; i < nX + 2; i++) {
    rowIB[i] = rowIB[i - 1] + 3;
  }
  rowIB[nX + 2] = rowIB[nX + 1];

  //          Matrix C
  j = 0;
  for (int i = 0; i < 3*nX; i+=3) {
    valC[i] = r2;             colC[i] = j++;
    valC[i + 1] = 1 - 2*r2;   colC[i + 1] = j++;
    valC[i + 2] = r2;         colC[i + 2] = j--;
  }

  rowIC[0] = 0;
  rowIC[1] = 1;
  for (int i = 2; i < nX + 2; i++) {
    rowIC[i] = rowIC[i - 1] + 3;
  }
  rowIC[nX + 2] = rowIC[nX + 1];


  // -----------------------------------------------------------------------
  //                         Вычисления
  //------------------------------------------------------------------------

  vector UNext(nX + 2);
  vector k1, k2, k3, k4;

  SpareMatrix A(valA, colA, rowIA, nX * 3, nX + 2);
  SpareMatrix B(valB, colB, rowIB, nX * 3, nX + 2);
  SpareMatrix C(valC, colC, rowIC, nX * 3, nX + 2);
  double g = dt/6.0;

  printf("%d\n\n", rowIB[nX+2]);

  double t0 = omp_get_wtime();
  for (int i =1; i <= sizeTime; i++) {
    k1 = A * U;
    k2 = B * k1;
    k3 = B * k2;
    k4 = C * k3;

    UNext = U + (k1 + k2*2 + k3*2 + k4)*g;
    std::swap (UNext, U);
  }
  double t1 = omp_get_wtime();

  //------------------------------------------------------------------------
  //                       Вывод результатов
  //------------------------------------------------------------------------

  printf("finish!\n");
  printf("time - %.15lf \n", t1 - t0);
  final(nX, U);

  //------------------------------------------------------------------------
  //                        Чистка памяти
  //------------------------------------------------------------------------

  delete[] valA;
  delete[] colA;
  delete[] rowIA;

  delete[] valB;
  delete[] colB;
  delete[] rowIB;

  return 0;

}

/*
____________________________________________________________________________


                          РЕАЛИЗАЦИЯ ФУНКЦИЙ

 ____________________________________________________________________________

*/


int init(double& xStart, double& xEnd, double& sigma, size_t& nX,
         double& tStart, double& tFinal, double& dt, int& check, vector& U) {
  FILE *fp;
  if ((fp = fopen("./../../../../initial/INPUT.txt", "r")) == NULL) {
    printf("Не могу найти файл!\n");
    return -1;
  }
  fscanf(fp, "XSTART=%lf\n", &xStart);
  fscanf(fp, "XEND=%lf\n", &xEnd);
  fscanf(fp, "SIGMA=%lf\n", &sigma);
  fscanf(fp, "NX=%lu\n", &nX);
  fscanf(fp, "TSTART=%lf\n", &tStart);
  fscanf(fp, "TFINISH=%lf\n", &tFinal);
  fscanf(fp, "dt=%lf\n", &dt);
  fscanf(fp, "BC=%d\n", &check);
  U.reserve(nX + 2);
  U.resize(nX + 2);

  // Заполнение функции в нулевой момент времени
  for (int i = 1; i < nX - 1; i++)
    fscanf(fp, "%lf", &U[i]);
  U[0] = U[nX + 1] = 0.0;
  fclose(fp);
  return 0;
}

int final(const size_t nX, vector& UFin) {
  FILE *fp;
  fp = fopen("./../../../../result/kirillRungeKuttaSparse.txt", "w");

  for (int i = 1; i < nX + 1; i++)
    fprintf(fp, "%.15le\n", UFin[i]);

  fclose(fp);
}
