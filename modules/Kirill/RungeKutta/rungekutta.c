#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main() {
  FILE *fp;

  int sizeTime = 0;
  int curTime, prevTime;

  const int N = 2;
  double* U[N];

  double xStart = 0.0, xEnd = 0.0;
  double sigma = 0.0;
  int nX = 1;
  double tStart = 0.0, tFinal = 0.0;
  double dt = 0.0;
  unsigned char check = 0;

  //------------------------------------------------------------------------
  //                Инициализация из файла
  //------------------------------------------------------------------------

  if((fp=fopen("./../../../../initial/INPUT.txt", "r")) == NULL) {
    printf("Не могу найти файл!\n");
    exit(-1);
  }

  // Заполнение конфигурационных настроек
  fscanf(fp,"XSTART=%lf\n", &xStart);
  fscanf(fp,"XEND=%lf\n", &xEnd);
  fscanf(fp,"SIGMA=%lf\n", &sigma);
  fscanf(fp,"NX=%d\n", &nX); // Не читает!?
  fscanf(fp,"TSTART=%lf\n", &tStart);
  fscanf(fp,"TFINISH=%lf\n", &tFinal);
  fscanf(fp,"dt=%lf\n", &dt);
  fscanf(fp,"BC=%hhu\n", &check);

  sizeTime = (int)((tFinal - tStart)/dt);

  printf("TIMESIZE = %d; NX = %d\n", sizeTime, nX);

  for (int x = 0; x < N; x++) {
    U[x] = (double *)malloc((nX + 2)*sizeof(double));
  }

  //  Вычисление шага по x
  double step = fabs(xStart - xEnd)/nX;

  if (2*sigma*dt > step*step) {
    printf("Выбор шага по времени не возможен в силу условия устойчивости!");
    exit(-1);
  }

  // Заполнение функции в нулевой момент времени
  for (int x = 1; x < nX -1; x++) {
    fscanf(fp, "%lf", &U[0][x]);
  }

  fclose(fp);

  // Задание граничных условий
  if (check == 2) {
    U[0][0] = U[0][1];
    U[0][nX - 1] = U[0][nX -2];
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

  double t0 = omp_get_wtime();
  for (int t = 1; t <= sizeTime; t++) {
    curTime = t%N;
    prevTime = (t + 1)%N;

    // Вычисление ki в данный момент времени
    for (int x = 1; x < nX + 2; x++)
      k1[x] = (U[prevTime][x + 1] - 2.0 * U[prevTime][x] + U[prevTime][x - 1])*coeff;

    for (int x = 1; x < nX + 2; x++)
      k2[x] = (U[prevTime][x + 1] + k1[x + 1]*dt*0.5 - 2.0*U[prevTime][x] -
          k1[x] * dt + U[prevTime][x - 1] + k1[x - 1] * dt * 0.5)*coeff;

    for (int x = 1; x < nX + 2; x++)
      k3[x] = (U[prevTime][x + 1] + k2[x + 1]*dt*0.5 - 2.0*U[prevTime][x] -
          k2[x]*dt + U[prevTime][x - 1] + k2[x - 1]*dt*0.5)*coeff;


    for (int x = 1; x < nX + 2; x++)
      k4[x] = (U[prevTime][x + 1] + k3[x + 1]*dt - 2.0*U[prevTime][x] -
          k3[x]*dt*2.0 + U[prevTime][x - 1] + k3[x - 1]*dt)*coeff;

    for (int x = 1; x < nX + 2; x++)
      U[curTime][x] = U[prevTime][x] + (k1[x] + 2.0*k2[x] + 2.0*k3[x] + k4[x])*h;

    // Задание граничных условий
    if (check == 2) {
      U[curTime][0] = U[curTime][1];
      U[curTime][nX - 1] = U[curTime][nX -2];
    } else if (check == 1) {
      // ???
      printf("HZ");
    }

  }

  double t1 = omp_get_wtime();
  printf("finish!\n");

  //------------------------------------------------------------------------
  //                       Вывод результатов и чистка памяти
  //------------------------------------------------------------------------

  printf("time - %.15lf \n", t1 - t0);
  fp = fopen("./../../../../result/kirillRungeKutta.txt", "w");
  for (int x = 1; x < nX + 1; x++)
    fprintf(fp, "%.15le\n", U[sizeTime%2][x]);
  fclose(fp);
  for (int i = 0; i < 2; i++)
    free(U[i]);

  return 0;
}