#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <complex.h>

// Complex по умолчанию double!
// OpenMP 4.0
#define ENABLE_PARALLEL 1

int main() {
  FILE *fp;

  int sizeTime = 0;
  int curTime, prevTime;

  const int N = 2;
  complex *U[N];

  complex xStart, xEnd;
  double sigma = 0.0;
  int nX = 1;
  double tStart = 0.0, tFinal = 0.0;
  double dt = 0.0;
  unsigned char check = 0;
  double step;

  //------------------------------------------------------------------------
  //                  ИНИЦИАЛИЗАЦИЯ И ЧТЕНИЕ ДАННЫХ
  //------------------------------------------------------------------------

  if ((fp = fopen("./../../../../initial/INPUT_COMPLEX.txt", "r")) == NULL) {
    printf("Не могу найти файл!\n");
    exit(-1);
  }

  double re, im = 0.;

  // Заполнение конфигурационных настроек
  fscanf(fp, "XSTART=%lf%lfi\n", &re, &im);
  xStart = CMPLX(re, im);

  fscanf(fp, "XEND=%lf%lfi\n", &re, &im);
  xEnd = CMPLX(re, im);

  fscanf(fp, "SIGMA=%lf\n", &sigma);
  fscanf(fp, "NX=%d\n", &nX);
  fscanf(fp, "TSTART=%lf\n", &tStart);
  fscanf(fp, "TFINISH=%lf\n", &tFinal);
  fscanf(fp, "dt=%lf\n", &dt);
  fscanf(fp, "BC=%hhu\n", &check);

  sizeTime = (int) ((tFinal - tStart) / dt);

  printf("XSTART=%lf%+lfi; XEND=%lf%+lfi; SIGMA=%lf; "
             "NX=%d; TSTART=%lf; TFINISH=%lf;"" dt=%lf; BC=%d;\n",
         creal(xStart),
         cimag(xStart),
         creal(xEnd),
         cimag(xEnd),
         sigma,
         nX,
         tStart,
         tFinal,
         dt,
         check);
  printf("TIMESIZE = %d; NX = %d\n", sizeTime, nX);

  for (int i = 0; i < N; i++) {
    U[i] = (complex *) malloc((nX + 2) * sizeof(complex));
  }

  //  Вычисление шага по x
  step = fabs(creal(xStart) - creal(xEnd)) / nX;
  if (2 * sigma * dt > step * step) {
    printf("Выбор шага по времени не возможен в силу условия устойчивости!\n");
    printf("%.10lf > %.10lf\n", 2 * sigma * dt, step * step);
    printf("Предлагаю взять dt = %.10lf\n", step * step / (2.0 * sigma));

    exit(-1);
  }

  // Заполнение функции в нулевой момент времени
  for (int i = 1; i < nX - 1; i++) {
    fscanf(fp, "%lf%lfi", &re, &im);
    U[0][i] = CMPLX(re, im);
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
  //                           ВЫЧИСЛЕНИЯ
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
    }
  double t1 = omp_get_wtime();

  //------------------------------------------------------------------------
  //                    ВЫВОД РЕЗУЛЬТАТОВ И ЧИСТКА МУСОРА
  //------------------------------------------------------------------------

  double diffTime = t1-t0;
  printf("finish!\n");
  printf("time run %.15lf\n", diffTime);
  printf("Gflops %.15lf\n", 5.0*nX*sizeTime/(diffTime*1000000000));

  fp = fopen("./../../../../result/kirillEulerComlplex.txt", "w");

  // Вывод результатов
  for (int i = 1; i < nX + 1; i++)
    fprintf(fp, "%.15lf%+.15lfi\n", creal(U[sizeTime%2][i]), cimag(U[sizeTime%2][i]));

  fclose(fp);

  // Чистка мусора
  for (int i = 0; i < 2; i++) {
    free(U[i]);
  }

  return 0;
}
