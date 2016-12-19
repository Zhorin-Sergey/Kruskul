#ifndef SOLUTIONS_KRUSKUL_SET_SET_H_
#define SOLUTIONS_KRUSKUL_SET_SET_H_

#include "stdio.h"

class Sepset {
 public:
  int *val;
  int size;
  Sepset() {}
  explicit Sepset(int size_);
  void singleton(int i);
  void merge(int i, int j);
  int poisk(int i);
  void clear();
  ~Sepset();
};

#endif  // SOLUTIONS_KRUSKUL_SET_SET_H_
