#include "Graph.h"

void reb(int a, int b, int s, Rebro *r) {
  r->a = a;
  r->b = b;
  r->size = s;
}

void gr(int vn, int rn, Rebro *rebr, Graph *g) {
  g->rebnumb = rn;
  g->vernumb = vn;
  g->rebra = new Rebro[rn];
  for (int i = 0; i < rn; i++)
    g->rebra[i] = rebr[i];
  for (int i = 0; i < rn; i++)
    for (int j = rn - 1; j > i; j--)
      if (g->rebra[j - 1].size < g->rebra[j].size) {
        Rebro tmp = g->rebra[j - 1];
        g->rebra[j - 1] = g->rebra[j];
        g->rebra[j] = tmp;
      }
}

void printt(Graph *g) {
  for (int i = 0; i < g->rebnumb; i++)
    std::cout << g->rebra[i].a << " - " << g->rebra[i].b << "; ";
}

void Krusk(Graph *g, Graph *ostov) {
  std::stack <Rebro> st;
  Sepset sss(g->vernumb);
  int n = g->vernumb;
  ostov->vernumb = n;
  ostov->rebnumb = n - 1;
  ostov->rebra = new Rebro[n - 1];
  for (int i = 0; i < n; i++)
    sss.singleton(i);
  for (int i = 0; i < g->rebnumb; i++)
    st.push(g->rebra[i]);
  int t = 0;
  while ((st.empty() != 1) && t < n - 1) {
    Rebro tmp = st.top();
    st.pop();
    if (sss.poisk(tmp.a) != sss.poisk(tmp.b)) {
      ostov->rebra[t] = tmp;
      t++;
      sss.merge(tmp.a, tmp.b);
    }
  }
}
