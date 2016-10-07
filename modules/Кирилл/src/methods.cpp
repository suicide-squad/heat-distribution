//
// Created by kirill on 07.10.16.
//

#include "methods.h"

#include "stdio.h"
#include "math.h"

#include <stdexcept>
#include <algorithm>


Methods::Methods(char* filename) {
  FILE* fp;
  if((fp=fopen(filename, "r")) == NULL)
    throw std::runtime_error("не могу найти файл");

  // Заполнение конфигурационных настроек
  fscanf(fp,"XSTART=%lf\n", &xStart);
  fscanf(fp,"XEND=%lf\n", &xEnd);
  fscanf(fp,"SIGMA=%lf\n", &sigma);
  fscanf(fp,"NX=%d\n", &nX); // Не читает!?
  fscanf(fp,"TSTART=%lf\n", &tStart);
  fscanf(fp,"TFINISH=%lf\n", &tFinal);
  fscanf(fp,"dt=%lf\n", &dt);
  fscanf(fp,"BC=%hhu\n", &check);

  sizeTime = static_cast<int>((tFinal - tStart)/dt);
  printf("%lf; %lf; %lf; %d; %lf; %lf; %lf; %d;\n",xStart, xEnd, sigma, nX, tStart, tFinal, dt, check);

  for (int i = 0; i < N; i++)
    U[i] = new double[nX + 2];

  startU = new double[nX];

  //  Вычисление шага по x
  step = fabs(xStart - xEnd)/nX;

  if (2*sigma*dt > step*step)
    throw std::runtime_error("выбор шага по времени не возможен в силу условия устойчивости");

  // Заполнение функции в нулевой момент времени
  for (int i = 0; i < nX; i++)
    fscanf(fp, "%lf", &startU[i]);


  fclose(fp);
}

Methods::~Methods() {
  for (int i = 0; i < N; i++)
    delete[] U[i];
  delete[] startU;
}

void Methods::run() {

  std::copy(startU, startU + nX, U[0] + 1);

  // Задание граничных условий
  if (check == 2) {
    U[0][0] = U[0][1];
    U[0][nX - 1] = U[0][nX -2];
  } else if (check == 1)
    // ???
    printf("HZ");

    int curTime, prevTime;

    // Заполнение сетки
    for (int i = 1; i <= sizeTime; i++) {
      curTime = i % N;
      prevTime = (i + 1) % N;

      for (int j = 1; j < nX + 2; j++)
        U[curTime][j] = (sigma * dt / (step * step)) * (U[prevTime][j + 1] -
            2 * U[prevTime][j] + U[prevTime][j - 1]) + U[prevTime][j];

      // Задание граничных условий
      if (check == 2) {
        U[curTime][0] = U[curTime][1];
        U[curTime][nX - 1] = U[curTime][nX - 2];
      }
    }
}

void Methods::saveResult(char* filename) {
  // Вывод результатов
  FILE* fp;

  fp = fopen(filename,"w");
  for (int i = 1; i < nX + 1; i++)
    fprintf(fp, "%2.15le\n", U[sizeTime%2][i]);

  fclose(fp);
}