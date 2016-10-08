//
// Created by kirill on 08.10.16.
//

#ifndef HEATDISTRIBUTION_EULER_H
#define HEATDISTRIBUTION_EULER_H

#include "methods.h"

class Euler: public Methods {
 public:
  Euler(char* filename):Methods(filename) {}
 protected:
  virtual double calculation(int t, int x) override;
};

#endif //HEATDISTRIBUTION_EULER_H
