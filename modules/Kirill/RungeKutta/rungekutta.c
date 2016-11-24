#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define ENABLE_PARALLEL 1

int main() {
  FILE *fp;


  const int N = 2;
  double* U;
  double* Unext;

  double xStart = 0.0, xEnd = 0.0;
  double sigma = 0.0;

  size_t sizeTime, nX;

  double tStart = 0.0, tFinal = 0.0;
  double dt = 0.0;
  unsigned char check = 0;

  //------------------------------------------------------------------------
  //                Инициализация из файла
  //------------------------------------------------------------------------

  if((fp=fopen("./../../initial/INPUT.txt", "r")) == NULL) {
    printf("Не могу найти файл!\n");
    exit(-1);
  }

  // Заполнение конфигурационных настроек
  fscanf(fp,"XSTART=%lf\n", &xStart);
  fscanf(fp,"XEND=%lf\n", &xEnd);
  fscanf(fp,"SIGMA=%lf\n", &sigma);
  fscanf(fp,"NX=%lu\n", &nX); // Не читает!?
  fscanf(fp,"TSTART=%lf\n", &tStart);
  fscanf(fp,"TFINISH=%lf\n", &tFinal);
  fscanf(fp,"dt=%lf\n", &dt);
  fscanf(fp,"BC=%hhu\n", &check);

  sizeTime = (size_t)((tFinal - tStart)/dt);

  #if ENABLE_PARALLEL
    printf("ПАРАЛЛЕЛЬНАЯ ВЕРСИЯ!\n");
  #endif
  printf("TIMESIZE = %lu; NX = %lu\n", sizeTime, nX);

  U = (double *)malloc((nX + 2)*sizeof(double));
  Unext = (double *)malloc((nX + 2)*sizeof(double));

  //  Вычисление шага по x
  double step = fabs(xStart - xEnd)/nX;

  if (2*sigma*dt > step*step) {
    printf("Выбор шага по времени не возможен в силу условия устойчивости!");
    exit(-1);
  }

  // Заполнение функции в нулевой момент времени
  for (int x = 1; x < nX -1; x++) {
    fscanf(fp, "%lf", &U[x]);
  }

  fclose(fp);

  // Задание граничных условий
  if (check == 2) {
    U[0] = U[1];
    U[nX + 1] = U[nX];
  } else if (check == 1) {
    // ???
    printf("HZ");
  }

  double* k1 =  (double *)malloc((nX + 2)*sizeof(double));
  double* k2 =  (double *)malloc((nX + 2)*sizeof(double));
  double* k3 =  (double *)malloc((nX + 2)*sizeof(double));
  double* k4 =  (double *)malloc((nX + 2)*sizeof(double));

  double coeff = 1.0/(step * step);
  double h = dt/6.0;

  // -----------------------------------------------------------------------
  //                         Вычисления
  //------------------------------------------------------------------------

  double *tmp;

  double t0 = omp_get_wtime();
  for (int t = 1; t <= sizeTime; t++) {
    // Вычисление ki в данный момент времени

    #pragma omp parallel for num_threads(2) if (ENABLE_PARALLEL)
    for (int x = 1; x < nX + 2; x++)
      k1[x] = (U[x + 1] - 2.0*U[x] + U[x - 1])*coeff;

    k1[0] = k1[1];
    k1[nX + 1] = k1[nX];

    #pragma omp parallel for num_threads(2) if (ENABLE_PARALLEL)
    for (int x = 1; x < nX + 2; x++)
      k2[x] = (U[x + 1] + k1[x + 1]*dt*0.5 - 2.0*U[x] -
          k1[x]*dt + U[x - 1] + k1[x - 1]*dt*0.5)*coeff;

    k2[0] = k2[1];
    k2[nX + 1] = k2[nX];

    #pragma omp parallel for num_threads(2) if (ENABLE_PARALLEL)
    for (int x = 1; x < nX + 2; x++)
      k3[x] = (U[x + 1] + k2[x + 1]*dt*0.5 - 2.0*U[x] -
          k2[x]*dt + U[x - 1] + k2[x - 1]*dt*0.5)*coeff;

    k3[0] = k3[1];
    k3[nX + 1] = k3[nX];

    #pragma omp parallel for num_threads(2) if (ENABLE_PARALLEL)
    for (int x = 1; x < nX + 2; x++)
      k4[x] = (U[x + 1] + k3[x + 1]*dt - 2.0*U[x] -
          k3[x]*dt*2.0 + U[x - 1] + k3[x - 1]*dt)*coeff;

    k4[0] = k4[1];
    k4[nX + 1] = k4[nX];

    #pragma omp parallel for num_threads(2) if (ENABLE_PARALLEL)
    for (int x = 1; x < nX + 2; x++)
      Unext[x] = U[x] + (k1[x] + 2.0*k2[x] + 2.0*k3[x] + k4[x])*h;

    // Задание граничных условий
    if (check == 2) {
      Unext[0] = Unext[1];
      Unext[nX + 1] = Unext[nX];
    } else if (check == 1) {
      // ???
      printf("HZ");
    }

    tmp = U;
    U = Unext;
    Unext = tmp;
  }

  double t1 = omp_get_wtime();
  printf("\nfinish!\n\n");

  //------------------------------------------------------------------------
  //                       Вывод результатов и чистка памяти
  //------------------------------------------------------------------------

  double diffTime = t1 - t0;
  double gflop = (nX*4 + nX*12 + nX*12 + nX*11 + nX*7)*sizeTime*1.0/1000000000.0;
  printf("Time\t%.15lf\n", diffTime);
  printf("GFlop\t%.lf\n", gflop);
  printf("GFlop's\t%.15lf\n", gflop*1.0/diffTime);

//  fp = fopen("./../../../../result/kirillRungeKutta.txt", "w");
//  for (int x = 1; x < nX + 1; x++)
//    fprintf(fp, "%.15le\n", U[x]);
//  fclose(fp);
  free(U);
  free(Unext);

  free(k1);
  free(k2);
  free(k3);
  free(k4);

  return 0;
}