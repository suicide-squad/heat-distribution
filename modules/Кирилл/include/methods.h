//
// Created by kirill on 07.10.16.
//

#ifndef HEATDISTRIBUTION_METHODS_H
#define HEATDISTRIBUTION_METHODS_H

class Methods {
 public:
  Methods(char*);

  virtual ~Methods();

  virtual double calculation(int t, int x) = 0;
  void run();
  void saveResult(char*);

 protected:
  static const int N = 2;

  double* U[N];
  double* startU;

  double xStart, xEnd;
  double sigma;
  int nX;
  double tStart, tFinal;
  double dt;
  unsigned char check;

  double step;

  int sizeTime;

};

#endif //HEATDISTRIBUTION_METHODS_H
