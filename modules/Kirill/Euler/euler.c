#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define ENABLE_PARALLEL 0

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
  double step;

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

  printf("XSTART=%lf; XEND=%lf; SIGMA=%lf; NX=%d; TSTART=%lf;"" TFINISH=%lf;"" dt=%lf; BC=%d;\n",
         xStart, xEnd, sigma, nX, tStart, tFinal, dt, check);

  for (int i = 0; i < N; i++) {
    U[i] = (double *)malloc((nX + 2)*sizeof(double));
  }

  //  Вычисление шага по x
  step = fabs(xStart - xEnd)/nX;

  if (2*sigma*dt > step*step) {
    printf("Выбор шага по времени не возможен в силу условия устойчивости!");
    exit(-1);
  }

  // Заполнение функции в нулевой момент времени
  for (int i = 1; i < nX -1; i++) {
    fscanf(fp, "%lf", &U[0][i]);
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

  // Заполнение сетки
  // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ

  double factor = sigma*dt/(step*step);

  double t0 = omp_get_wtime();
  #pragma omp parallel if (ENABLE_PARALLEL)
  {
    for (int i = 1; i <= sizeTime; i++) {
      curTime = i % N;
      prevTime = (i + 1) % N;

      int j;
      #pragma omp for
      for (j = 1; j < nX + 2; j++)
        U[curTime][j] = factor * (U[prevTime][j + 1] - 2.0 * U[prevTime][j] +
            U[prevTime][j - 1]) + U[prevTime][j];

      // Задание граничных условий
      if (check == 2) {
        U[curTime][0] = U[curTime][1];
        U[curTime][nX - 1] = U[curTime][nX - 2];
      } else if (check == 1) {
        printf("HZ");
      }
    }
  }
  double t1 = omp_get_wtime();

  printf("finish!\n");
  printf("time run %.15lf\n", t1-t0);

  fp = fopen("./../../../../result/kirillEuler.txt", "w");

  // Вывод результатов
  for (int j = 1; j < nX + 1; j++)
    fprintf(fp, "%.15le\n", U[sizeTime%2][j]);

  fclose(fp);

  // Чистка мусора
  for (int i = 0; i < 2; i++) {
    free(U[i]);
  }

  return 0;
}