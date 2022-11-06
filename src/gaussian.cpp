#include "internal.hpp"
#include "util.hpp"
#include <immintrin.h>
#include <initializer_list>
#include <math.h>

namespace CaDiCaL {

/*----------------------------------------------------------------------*/
// Raw_XOR_Equations for equation partition and dispatch
Raw_XOR_Equations::Raw_XOR_Equations(Internal *internal) : internal(internal) {}

Raw_XOR_Equations::~Raw_XOR_Equations() {
  for (auto *cl : equations)
    delete cl;
}

/*----------------------------------------------------------------------*/

void Raw_XOR_Equations::add_clause(const vector<int> &lits, int flip) {
  Raw_XOR_Clause *e = new Raw_XOR_Clause;
  e->vars = lits;
  e->flip = flip;

  equations.push_back(e);
}

/*----------------------------------------------------------------------*/

int Raw_XOR_Equations::find(int x) {
  return parent[x] == x ? x : parent[x] = find(parent[x]);
}

/*----------------------------------------------------------------------*/

void Raw_XOR_Equations::merge(int x, int y) {
  int u = find(x), v = find(y);
  if (u == v)
    return;
  parent[u] = v;
  component_size[v] += component_size[u];
}

/*----------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------*/

void Raw_XOR_Equations::analyze_connectivity() {
  int eq_cnt = equations.size();
  parent.resize(eq_cnt, -1);
  component_size.resize(eq_cnt, 1);

  for (int i = 0; i < eq_cnt; ++i)
    parent[i] = i;
  for (int i = 0; i < eq_cnt; ++i)
    for (int j = i + 1; j < eq_cnt; ++j)
      if (check_incidence(equations[i]->vars, equations[j]->vars))
        merge(i, j);

  int component_cnt = 0;
  for (int i = 0; i < eq_cnt; ++i)
    if (find(i) == i)
      printf("Component %d with size %d\n", component_cnt++, component_size[i]);
}

/*----------------------------------------------------------------------*/

void Raw_XOR_Equations::partition() {
  // assume analyze_connectivity is has been run
  int eq_cnt = equations.size();
  int block_cnt = 0;

  vector<int> block_id(eq_cnt);

  for (int i = 0; i < eq_cnt; ++i)
    if (find(i) == i)
      block_id[i] = block_cnt++;

  vector<vector<Raw_XOR_Clause *>> equation_blocks(block_cnt);
  for (int i = 0; i < eq_cnt; ++i) {
    int blk = block_id[find(i)];
    equation_blocks[blk].push_back(equations[i]);
  }

  blocks.resize(block_cnt);
  for (int i = 0; i < block_cnt; ++i) {
    blocks[i] = new GaussianDense(internal, this, i);
    blocks[i]->initialize(equation_blocks[i]);
  }
}

/*----------------------------------------------------------------------*/

void Raw_XOR_Equations::analyze() {
  int size = equations.size();
  int var_cnt = 0;

  for (auto *e : equations)
    var_cnt += e->vars.size();

  printf("Num Lits: %d, Num Eqs: %d\n", var_cnt, size);

  analyze_connectivity();
  partition();
}

/*----------------------------------------------------------------------*/

vector<pair<int, int>> Raw_XOR_Equations::propagate(const vector<int> &lit,
                                                    int level) {
  START(xpropagate);
  int num_block = blocks.size();
  vector<vector<int>> var_by_block(num_block);
  vector<vector<int>> assgn_by_block(num_block);

  for (int i = 0, size = lit.size(); i < size; ++i) {
    int var = abs(lit[i]);
    if (internal->get_gaussian_block(var) != -1) {
      unsigned int assgn = (bign(lit[i]) >> 1) ^ 1; // bign(false)=2, bign(true)=1
      int blk = internal->get_gaussian_block(var);
      var_by_block[blk].push_back(var);
      assgn_by_block[blk].push_back(assgn);
    }
  }

  vector<pair<int, int>> deduced_lit;

  for (int i = 0; i < num_block; ++i)
    if (var_by_block[i].size() != 0) {
      auto vec =
	blocks[i]->propagate(var_by_block[i], assgn_by_block[i], level);
      deduced_lit.insert(deduced_lit.end(), vec.begin(), vec.end());
      if (conflict_flag) {
        generated_reasons.push_back(make_pair(conflict, level));
	STOP(xpropagate);
        return {};
      }
    }

  STOP(xpropagate);
  return deduced_lit;
}

/*----------------------------------------------------------------------*/

void Raw_XOR_Equations::backtrack(int level) {
  START(xbacktrack);
  LOG("enter xbacktrack level %d\n", level);
  conflict_flag = false;
  for (auto *blk : blocks)
    blk->backtrack(level);

  auto i = generated_reasons.begin();
  auto end = generated_reasons.end();
  for (auto j = i; j != end; ++j) {
    if (j->second > level)
      delete[](char*)(j->first);
    else
      *i++ = *j;
  }
  generated_reasons.resize(i - generated_reasons.begin());
  STOP(xbacktrack);
  LOG("exit xbacktrack");
}

Clause *Raw_XOR_Equations::lazily_gen_reason(int block, int eqn) {
  auto reason = blocks[block]->analyze(eqn);
  generated_reasons.push_back(reason);
  return reason.first;
}

bool Raw_XOR_Equations::all_satisfied() {
  for (int i = 0; i < blocks.size(); ++i)
    if (!blocks[i]->all_satisfied())
      return false;
  return true;
}

/*----------------------------------------------------------------------*/

}; // namespace CaDiCaL
