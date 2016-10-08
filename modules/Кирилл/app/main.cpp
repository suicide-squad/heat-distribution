#include "euler.h"

int main() {
  Euler m("./../../../initial/INPUT.txt");
  m.run();
  m.saveResult("./../../../result/PetrovResult.txt");

  return 0;
}