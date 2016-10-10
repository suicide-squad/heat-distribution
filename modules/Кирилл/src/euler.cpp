//
// Created by kirill on 08.10.16.
//

#include "euler.h"

double Euler::calculation(int t, int x) {
  return (sigma*dt/(step*step))*(U[t][x + 1] -
      2*U[t][x] + U[t][x - 1]) + U[t][x];
}