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

  printf("%d",sizeTime);


  //  чистка и вывод результатов
  final(nX, U);
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
