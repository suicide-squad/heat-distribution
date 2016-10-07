#include "methods.h"

int main() {
  Methods m("./../../../initial/INPUT.txt");
  m.run();
  m.saveResult("./../../../result/PetrovResult.txt");

  return 0;
}