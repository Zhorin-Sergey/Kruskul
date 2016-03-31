#ifndef SOLUTIONS_KRUSKUL_SET_GRAPH_H_
#define SOLUTIONS_KRUSKUL_SET_GRAPH_H_

#include <iostream>
#include <stack>
#include "Set.h"

struct Rebro {
  int a;
  int b;
  int size;
};

struct Graph {
  int vernumb;
  int rebnumb;
  Rebro* rebra;
};

void reb(int, int, int, Rebro*);
void gr(int, int, Rebro*, Graph*);
void printt(Graph *g);
void Krusk(Graph*, Graph*);

#endif  // SOLUTIONS_KRUSKUL_SET_GRAPH_H_
