//
// Created by kirill on 24.10.16.
//

#include <stdio.h>
#include <math.h>

#include "sp_mat.h"

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
  // инициализация данных
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

  double h = sigma*dt/(step*step);
  vector UNext(nX + 2);

  double* val = new double[nX*3 + 2];
  int* col = new int[nX*3 + 2];
  int* rowI = new int[nX + 2 + 1];

  //------------------------------------------------------------------------
  // Заполнение значений и номера столбцов матрицы
  //------------------------------------------------------------------------

  val[0] = 1; col[0] = 0;
  int j = 0;
  for (int i = 1; i < 3*nX + 1; i+=3) {
    val[i] = h; col[i] = j++;
    val[i + 1] = 1 - 2*h; col[i + 1] = j++;
    val[i + 2] = h; col[i + 2] = j--;
  }
  val[nX*3 + 1] = 1; col[3*nX + 1] = static_cast<int>(nX) + 1;

  rowI[0] = 0;
  rowI[1] = 1;
  for (int i = 2; i < nX + 2; i++)
    rowI[i] = rowI[i - 1] + 3;
  rowI[nX + 2] = rowI[nX + 1] + 1;

  // -----------------------------------------------------------------------
  // Вычисления
  //------------------------------------------------------------------------

  SpareMatrix A(val, col, rowI, nX * 3 + 2, nX + 2);

  double t0 = omp_get_wtime();
  for (int i =1; i <= sizeTime; i++) {
    UNext = A * U;
    std::swap (UNext, U);
  }
  double t1 = omp_get_wtime();

  //------------------------------------------------------------------------
  // Вывод результатов
  //------------------------------------------------------------------------

  printf("finish!\n");
  printf("time - %.15lf \n", t1 - t0);
  final(nX, U);
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
  fp = fopen("./../../../../result/kirillEulerSparse.txt", "w");

  for (int i = 1; i < nX + 1; i++)
    fprintf(fp, "%.15le\n", UFin[i]);

  fclose(fp);
}
