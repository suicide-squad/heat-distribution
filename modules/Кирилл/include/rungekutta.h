//
// Created by kirill on 08.10.16.
//

#ifndef HEATDISTRIBUTION_RUNGEKUTTA_H
#define HEATDISTRIBUTION_RUNGEKUTTA_H

#include "methods.h"

class RungeKutta : public Methods {
 public:
  RungeKutta(char* filename):Methods(filename) {}
 protected:
  virtual double calculation(int t, int x) override;

};

#endif //HEATDISTRIBUTION_RUNGEKUTTA_H
