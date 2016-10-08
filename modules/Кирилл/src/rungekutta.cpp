//
// Created by kirill on 08.10.16.
//

#include "rungekutta.h"

double RungeKutta::calculation(int t, int x) {
  double k1,k2,k3,k4;

  k1 = sigma*(U[t][x + 1] - 2*U[t][x] + U[t][x - 1])/(step*step);

  double U1 = U[x][t] + dt*k1/2;

 // k2 = ???

  return U[t][x] + step*(k1 + 2*k2 + 2*k3 + k4)/6;
}