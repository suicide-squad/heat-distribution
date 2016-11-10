#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// OpenMP 4.0
#define ENABLE_PARALLEL 0

int main() {
  FILE *fp;
  const int N = 2;

  int curTime, prevTime;

  double *U[N];

  size_t nX, sizeTime;

  double xStart = 0.0, xEnd = 0.0;
  double sigma = 0.0;
  double tStart = 0.0, tFinal = 0.0;
  double dt = 0.0;
  unsigned char check = 0;
  double step;

  //------------------------------------------------------------------------
  //                      ИНИЦИАЛИЗАЦИЯ И ЧТЕНИЕ ДАННЫХ
  //------------------------------------------------------------------------

  if ((fp = fopen("./../../../../initial/INPUT.txt", "r")) == NULL) {
    printf("Не могу найти файл!\n");
    exit(-1);
  }

  // Заполнение конфигурационных настроек
  fscanf(fp, "XSTART=%lf\n", &xStart);
  fscanf(fp, "XEND=%lf\n", &xEnd);
  fscanf(fp, "SIGMA=%lf\n", &sigma);
  fscanf(fp, "NX=%lu\n", &nX);
  fscanf(fp, "TSTART=%lf\n", &tStart);
  fscanf(fp, "TFINISH=%lf\n", &tFinal);
  fscanf(fp, "dt=%lf\n", &dt);
  fscanf(fp, "BC=%hhu\n", &check);

  sizeTime = (size_t) ((tFinal - tStart) / dt);

  printf("TIMESIZE = %lu; NX = %lu\n", sizeTime, nX);

  for (int i = 0; i < N; i++) {
    U[i] = (double *) malloc((nX + 2) * sizeof(double));
  }

  //  Вычисление шага по x
  step = fabs(xStart - xEnd) / nX;
  if (2 * sigma * dt > step * step) {
    printf("Выбор шага по времени не возможен в силу условия устойчивости!\n");
    printf("%.10lf > %.10lf\n", 2 * sigma * dt, step * step);
    printf("Предлагаю взять dt = %.10lf\n", step * step / (2.0 * sigma));

    exit(-1);
  }

  // Заполнение функции в нулевой момент времени
  for (int i = 1; i < nX - 1; i++) {
    fscanf(fp, "%lf", &U[0][i]);
  }

  fclose(fp);

  // Задание граничных условий
  if (check == 2) {
    U[0][0] = U[0][1];
    U[0][nX - 1] = U[0][nX - 2];
  } else if (check == 1) {
    // ???
    printf("HZ");
  }

  // Заполнение сетки
  //------------------------------------------------------------------------
  //                      ОСНОВНЫЕ ВЫЧИСЛЕНИЯ
  //------------------------------------------------------------------------

  double factor = sigma * dt / (step * step);

  int j;
  double t0 = omp_get_wtime();

  for (int i = 1; i <= sizeTime; i++) {
    curTime = i % N;
    prevTime = (i + 1) % N;


    #pragma parallel omp num_threads(2) if (ENABLE_PARALLEL)
    {
      #pragma omp for nowait
      for (j = 1; j < nX + 2; j++)
        U[curTime][j] = factor * (U[prevTime][j + 1] - 2.0 * U[prevTime][j] +
            U[prevTime][j - 1]) + U[prevTime][j];
    }
    // Задание граничных условий
    if (check == 2) {
      U[curTime][0] = U[curTime][1];
      U[curTime][nX - 1] = U[curTime][nX - 2];
    } else if (check == 1) {
      printf("HZ");
    }
  } // for
  double t1 = omp_get_wtime();
  printf("finish!\n\n");

  //------------------------------------------------------------------------
  //                      ВЫВОД РЕЗУЛЬТАТОВ И ЧИСТКА МУСОРА
  //------------------------------------------------------------------------

  double diffTime = t1-t0;
  unsigned long long flop = 5*nX*sizeTime;
  printf("Time\t%.15lf\n", diffTime);
  printf("Flop\t%.0llu\n", flop);
  printf("GFlops\t%.15lf\n", 1.0*flop/(diffTime*1000000000.0));



  fp = fopen("./../../../../result/kirillEuler.txt", "w");

  // Вывод результатов
  for (int i = 1; i < nX + 1; i++)
    fprintf(fp, "%.15le\n", U[sizeTime%2][i]);

  fclose(fp);

  // Чистка мусора
  for (int i = 0; i < 2; i++) {
    free(U[i]);
  }

  return 0;
}