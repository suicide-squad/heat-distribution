#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main() {
  FILE *fp;

  int i, j;
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
  double step = 0.0;

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

  nX = 500;
  sizeTime = (int)((tFinal - tStart)/dt);

  printf("%lf; %lf; %lf; %d; %lf; %lf; %lf; %d;\n",xStart, xEnd, sigma, nX, tStart, tFinal, dt, check);

  for (i = 0; i < N; i++) {
    U[i] = (double *)malloc((nX + 2)*sizeof(double));
  }

  //  Вычисление шага по x
  step = fabs(xStart - xEnd)/nX;

  if (2*sigma*dt > step*step) {
    printf("Выбор шага по времени не возможен в силу условия устойчивости!");
    exit(-1);
  }

  // Заполнение функции в нулевой момент времени
  for (i = 1; i < nX -1; i++) {
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

  clock_t t0 = clock();
  for (i = 1; i <= sizeTime; i++) {
    curTime = i%N;
    prevTime = (i + 1)%N;

    for (j = 1; j < nX + 2; j++)
      U[curTime][j] = factor*(U[prevTime][j + 1] - 2*U[prevTime][j] +
                      U[prevTime][j - 1]) + U[prevTime][j];

    // Задание граничных условий
    if (check == 2) {
      U[curTime][0] = U[curTime][1];
      U[curTime][nX - 1] = U[curTime][nX -2];
    } else if (check == 1) {
      printf("HZ");
    }

  }
  clock_t t1 = clock();

  printf("finish!\n");
  printf("time run %ld\n", t1-t0);

  fp = fopen("./../../../../result/kirillEuler.txt", "w");

  // Вывод результатов
  for (j = 1; j < nX + 1; j++)
    fprintf(fp, "%.15le\n", U[sizeTime%2][j]);

  fclose(fp);

  // Чистка мусора
  for (i = 0; i < 2; i++) {
    free(U[i]);
  }

  return 0;
}