#include "Graph.h"

void main() {
  Rebro a;
  reb(0, 1, 1, &a);
  Rebro b;
  reb(1, 2, 4, &b);
  Rebro c;
  reb(2, 3, 5, &c);
  Rebro d;
  reb(3, 4, 3, &d);
  Rebro e;
  reb(1, 3, 1, &e);
  Rebro f;
  reb(1, 4, 1, &f);
  Rebro g;
  reb(0, 4, 2, &g);
  Rebro rrr[7] = { a, b, c, d, e, f, g };
  Graph ggg;
  gr(5, 7, rrr, &ggg);
  Graph ooo;
  Krusk(&ggg, &ooo);
  printt(&ooo);
}
