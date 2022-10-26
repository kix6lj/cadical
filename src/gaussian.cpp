#include "internal.hpp"

namespace CaDiCaL {
Raw_XOR_Equations::Raw_XOR_Equations(Internal *internal) :
  internal(internal)
{}

Raw_XOR_Equations::~Raw_XOR_Equations() {
  for (auto *cl : equations)
    delete cl;
}

void Raw_XOR_Equations::add_clause(const vector<int> &lits, int flip) {
  Raw_XOR_Clause *e = new Raw_XOR_Clause;
  e->vars = lits;
  e->flip = flip;

  equations.push_back(e);
}

int Raw_XOR_Equations::find(int x) {
  return parent[x] == x ? x : parent[x] = find(parent[x]);
}

void Raw_XOR_Equations::merge(int x, int y) {
  int u = find(x), v = find(y);
  if (u == v)
    return;
  parent[u] = v;
  component_size[v] += component_size[u];
}

bool Raw_XOR_Equations::check_incidence(const vector<int> &vars1,
                                        const vector<int> &vars2) {
  bool flag = false;
  for (auto v : vars1)
    internal->mark(v);
  
  for (auto v : vars2)
    if (internal->marked(v)) {
      flag = true;
      break;
    }
  for (auto v : vars1)
    internal->unmark(v);
  return flag;
}

void Raw_XOR_Equations::analyze_connectivity() {
  size_t eq_cnt = equations.size();
  size_t var_cnt = 0;
  parent.resize(eq_cnt, -1);
  component_size.resize(eq_cnt, 1);

  for (size_t i = 0; i < eq_cnt; ++i)
    parent[i] = i;
  for (size_t i = 0; i < eq_cnt; ++i)
    for (size_t j = i + 1; j < eq_cnt; ++j)
      if (check_incidence(equations[i]->vars, equations[j]->vars))
        merge(i, j);

  int component_cnt = 0;
  for (int i = 0; i < eq_cnt; ++i)
    if (parent[i] == i)
      printf("Component %d with size %d\n", component_cnt++,
             component_size[i]);
}

void Raw_XOR_Equations::analyze() {
  size_t size = equations.size();
  size_t var_cnt = 0;

  for (auto *e : equations)
    var_cnt += e->vars.size();

  printf("Num Lits: %d, Num Eqs: %d\n", var_cnt, size);

  analyze_connectivity();
}
}; // namespace CaDiCaL
