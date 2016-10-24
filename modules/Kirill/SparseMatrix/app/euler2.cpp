//
// Created by kirill on 24.10.16.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sparse_matrix.h>

int init(double&, double&, double&, int&, double&, double&, double&, int&, double*&);
int final(const int, double*&);

int main() {
  double xStart, xEnd;
  double sigma;
  int nX;
  double tStart, tFinal;
  double dt;
  int check;

  double* U;

  // инициализация данных из файла
  init(xStart, xEnd, sigma, nX, tStart, tFinal, dt, check, U);

  double step = fabs(xStart - xEnd) / nX;
  int sizeTime = static_cast<int>((tFinal - tStart) / dt);
  if (2 * sigma * dt > step * step) {
    printf("Выбор шага по времени не возможен в силу условия устойчивости!\n");
    printf("%.10lf > %.10lf\n", 2 * sigma * dt, step * step);
    printf("Предлагаю взять dt = %.10lf\n", step * step / (2.0 * sigma));

    return -1;
  }

  double h = sigma*dt/(step*step);
  double* UNext = new double [nX + 2];

  double* val = new double[nX*3 + 2];
  int* col = new int[nX*3 + 2];
  int* rowI = new int[nX + 2 + 1];


  // Заполнение значений и номера столбцов матрицы
  val[0] = 1; col[0] = 0;
  int j = 0;
  for (int i = 1; i < 3*nX + 1; i+=3) {
    val[i] = h; col[i] = j++;
    val[i + 1] = 1 - 2*h; col[i + 1] = j++;
    val[i + 2] = h; col[i + 2] = j--;
  }
  val[nX*3 + 1] = 1; col[3*nX + 1] = nX + 1;

  rowI[0] = 0;
  rowI[1] = 1;
  for (int i = 2; i < nX + 2; i++)
    rowI[i] = rowI[i - 1] + 3;
  rowI[nX + 2] = rowI[nX + 1] + 1;

  // TEST
//  printf("h = %.4lf\n",h);
//  for (int i = 0; i < nX + 3; i++)
//    printf("%d\n", rowI[i]);

  SpareMatrix A(val, col, rowI, nX*3 + 2, nX + 2);
  UNext = A * U;

  //  чистка и вывод результатов
  final(nX, UNext);
  return 0;

}


int init(double& xStart, double& xEnd, double& sigma, int& nX,
         double& tStart, double& tFinal, double& dt, int& check, double*& U) {
  FILE *fp;
  if ((fp = fopen("./../../../../initial/INPUT.txt", "r")) == NULL) {
    printf("Не могу найти файл!\n");
    return -1;
  }
  fscanf(fp, "XSTART=%lf\n", &xStart);
  fscanf(fp, "XEND=%lf\n", &xEnd);
  fscanf(fp, "SIGMA=%lf\n", &sigma);
  fscanf(fp, "NX=%d\n", &nX);
  fscanf(fp, "TSTART=%lf\n", &tStart);
  fscanf(fp, "TFINISH=%lf\n", &tFinal);
  fscanf(fp, "dt=%lf\n", &dt);
  fscanf(fp, "BC=%d\n", &check);
  U = new double [nX + 2];

  // Заполнение функции в нулевой момент времени
  for (int i = 1; i < nX - 1; i++) {
    fscanf(fp, "%lf", U + i);
  }
  fclose(fp);
  return 0;
}


int final(const int nX, double*& UFin) {
  FILE *fp;
  fp = fopen("./../../../../result/kirillEuler2.txt", "w");

  for (int i = 1; i < nX + 1; i++)
    fprintf(fp, "%.15le\n", UFin[i]);

  fclose(fp);
  delete[] UFin;
}
