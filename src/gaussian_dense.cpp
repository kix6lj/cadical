#include "internal.hpp"
#include <cstdio>
#include <cstdlib>
#include <immintrin.h>
#include <initializer_list>
#include <math.h>

// =------------------------------=
// GaussianDense

#define TEMP_EQN 2
#define WIDTH 256

#define DIV_CEIL(x, y) ((x) + (y)-1) / (y)
#define DIV_FLOOR(x, y) (x) / (y)

namespace CaDiCaL {
GaussianDense::GaussianDense(Internal *internal, Raw_XOR_Equations *manager,
                             int index)
    : internal(internal), manager(manager), index(index) {}

GaussianDense::~GaussianDense() {
  // delete[] equations;
  // delete[] assignment;
  // delete[] marks;
  // delete[] base_mask;
  // delete[] watch_mask;
  // delete[] temp_mem;
  free(equations);
  free(assignment);
  free(marks);
  free(base_mask);
  free(watch_mask);
  free(temp_mem);
}

void GaussianDense::init_watch() {
  base_col.resize(num_eqns);
  watch_col.resize(num_eqns);
  base_row.resize(num_vars + 1);
  col_is_base.resize(num_vars + 1, 0);
  col_is_watch.resize(num_vars + 1, 0);
  watched_eqns.resize(num_vars + 1);

  for (int i = 0; i < num_eqns; ++i) {
    char *ptr = eqn_start_ptr(i);
    int base = msb_position(ptr);

    if (base == 0) {
      manager->conflict_flag = true;
      continue;
    }

    base_col[i] = base;
    base_row[base] = i;

    // elimination
    for (int j = 0; j < num_eqns; ++j)
      if (i != j) {
        char *ptr_j = eqn_start_ptr(j);
        if (get_bit(ptr_j, base))
          XOR(ptr_j, ptr);
      }
  }

  for (int i = 0; i < num_eqns; ++i) {
    char *ptr = eqn_start_ptr(i);
    int base = base_col[i];

    assert(base != 0);
    unset_bit(ptr, base);
    watch_col[i] = msb_position(ptr);
    set_bit(ptr, base);
    assert(get_bit(ptr, watch_col[i]));

    if (watch_col[i] == 0) {
      // FIXME: needs propagate directly, although most of the time this won't
      // happen
      assert(0 && "Propagate at level 0 is not implemented");
      continue;
    }
  }

  for (int i = 0; i < num_eqns; ++i) {
    set_bit(base_mask, base_col[i]);
    set_bit(watch_mask, watch_col[i]);
    set_watch(watch_col[i], i);
  }
}

void GaussianDense::initialize(const vector<Raw_XOR_Clause *> &original) {
  for (auto *c : original)
    for (auto x : c->vars)
      if (!internal->marked(x)) {
        assert(x <= internal->max_var);
        var_id.push_back(x);
        internal->mark(x);
      }

  for (auto x : var_id)
    internal->unmark(x);

  LOG("Block has %d equations and %d vars\n", original.size(), var_id.size());

  num_vars = var_id.size();
  num_eqns = original.size();
  num_bytes_per_eqn = DIV_CEIL(num_vars + 1, WIDTH) * (WIDTH >> 3);

  // equations = new char[num_bytes_per_eqn * original.size()];
  // assignment = new char[num_bytes_per_eqn];
  // marks = new char[num_bytes_per_eqn];
  // base_mask = new char[num_bytes_per_eqn];
  // watch_mask = new char[num_bytes_per_eqn];
  // temp_mem = new char[num_bytes_per_eqn * TEMP_EQN];

  equations = (char *)aligned_alloc(32, num_bytes_per_eqn * original.size());
  assignment = (char *)aligned_alloc(32, num_bytes_per_eqn);
  marks = (char *)aligned_alloc(32, num_bytes_per_eqn);
  base_mask = (char *)aligned_alloc(32, num_bytes_per_eqn);
  watch_mask = (char *)aligned_alloc(32, num_bytes_per_eqn);
  temp_mem = (char *)aligned_alloc(32, num_bytes_per_eqn * TEMP_EQN);

  xup_level.resize(num_eqns);

  memset(equations, 0, num_bytes_per_eqn * original.size());
  memset(assignment, 0, num_bytes_per_eqn);
  memset(marks, 0, num_bytes_per_eqn);
  memset(base_mask, 0, num_bytes_per_eqn);
  memset(watch_mask, 0, num_bytes_per_eqn);

  for (int i = 0; i < var_id.size(); ++i) {
    internal->set_gaussian_pos(var_id[i], i + 1);
    internal->set_gaussian_block(var_id[i], index);
  }

  for (int i = 0; i < original.size(); ++i) {
    char *eqn_ptr = eqn_start_ptr(i);
    if (original[i]->flip)
      set_bit(eqn_ptr, 0);
    for (int x : original[i]->vars)
      set_bit(eqn_ptr, internal->get_gaussian_pos(x));
  }

  init_watch();
}

void GaussianDense::set_watch(int col, int eqn) {
  watched_eqns[col].push_back(eqn);
  watch_col[eqn] = col;
  set_bit(watch_mask, col);
}

void GaussianDense::remove_watch(int col, int eqn) {
  auto &we = watched_eqns[col];
  const auto end = we.end();
  auto i = we.begin();

  for (auto j = i; j != end; j++) {
    const int &w = *i++ = *j;
    if (w == eqn)
      i--;
  }

  LOG("REMOVE %d %d\n", var_id[col - 1], eqn);
  assert(i + 1 == end);
  we.resize(i - we.begin());
  if (we.size() == 0)
    unset_bit(watch_mask, col);
  watch_col[eqn] = 0;
}

vector<pair<int, int>> GaussianDense::propagate(const vector<int> &vars,
                                                const vector<int> &assgns,
                                                int level) {
  vector<int> base_cols;
  vector<int> watch_cols;
  vector<int> update_watch_rows;

  /*---------------------------------------------------------------------*/

  for (int i = 0; i < vars.size(); ++i) {
    int var = vars[i];
    int assgn = assgns[i];

    int pos = internal->get_gaussian_pos(var);
    assert(var > 0);
    assert(!get_bit(marks, pos));
    assignment[pos >> 3] |= assgn << (pos & 7);
    marks[pos >> 3] |= 1 << (pos & 7);

    if (get_bit(base_mask, pos))
      base_cols.push_back(pos);

    if (get_bit(watch_mask, pos)) {
      update_watch_rows.insert(update_watch_rows.end(),
                               watched_eqns[pos].begin(),
                               watched_eqns[pos].end());
      for (auto eqn : watched_eqns[pos])
        watch_col[eqn] = 0;
    }
  }

  assigned.insert(assigned.end(), vars.begin(), vars.end());

  /*----------------------------------------------------------------------*/

  // Update bases
  if (base_cols.size() > 0) {
    // perform elimination on each col's related row
    // each col belongs to exactly one row
    for (int col : base_cols) {
      int eqn = base_row[col];
      char *ptr = eqn_start_ptr(eqn);
      int next_base;

      if (watch_col[eqn] && !get_bit(marks, watch_col[eqn])) {
        next_base = watch_col[eqn];
        assert(get_bit(ptr, watch_col[eqn]));
        assert(!get_bit(marks, watch_col[eqn]));
        // watch was used by base, need to call remove_watch
        remove_watch(watch_col[eqn], eqn);
        update_watch_rows.push_back(eqn);
      } else {
        // ALLOCATE MEM
        char *remained_ones = temp_mem;
        IAND(remained_ones, ptr, marks);
        if (count(remained_ones) == 0) {
          // ALLOCATE MEM
          char *evaluate = temp_mem + num_bytes_per_eqn;
          AND(evaluate, ptr, assignment);
          int xor_sum = count(evaluate) & 1;
          int required = !!get_bit(ptr, 0);
          if (xor_sum != required) {
            manager->conflict = conflict(eqn);
            manager->conflict_flag = true;
            return {};
          }
          continue;
        }
        watch_col[eqn] = 0;

        next_base = msb_position(remained_ones);
        update_watch_rows.push_back(eqn);
      }

      for (int j = 0; j < num_eqns; ++j)
        if (eqn != j) {
          char *ptr_j = eqn_start_ptr(j);
          if (get_bit(ptr_j, next_base)) {
            XOR(ptr_j, ptr);
            assert(get_bit(ptr, next_base));
            assert(!get_bit(ptr_j, next_base));
            if (watch_col[j] && !get_bit(ptr_j, watch_col[j])) {
              remove_watch(watch_col[j], j);
              update_watch_rows.push_back(j);
            }
          }
        }

      assert(!get_bit(base_mask, next_base));
      unset_bit(base_mask, base_col[eqn]);
      set_bit(base_mask, next_base);
      base_col[eqn] = next_base;
      base_row[next_base] = eqn;
      unset_bit(base_mask, col);
    }
  }

  /*-------------------------------------------------------------------------*/

  sort(update_watch_rows.begin(), update_watch_rows.end());
  update_watch_rows.resize(
      unique(update_watch_rows.begin(), update_watch_rows.end()) -
      update_watch_rows.begin());

  /*-------------------------------------------------------------------------*/

  // Update watches
  vector<pair<int, int>> propagated_lits;
  if (update_watch_rows.size() > 0) {
    for (auto eqn : update_watch_rows) {

      if (watch_col[eqn] && !get_bit(marks, watch_col[eqn]))
        continue;

      watch_col[eqn] = 0;

      char *ptr = eqn_start_ptr(eqn);
      char *remained_ones = temp_mem;

      IAND(remained_ones, ptr, marks);
      int num_ones = count(remained_ones);

      if (num_ones == 1) {
        char *evaluate = temp_mem + num_bytes_per_eqn;
        AND(evaluate, ptr, assignment);
        int xor_sum = count(evaluate) & 1;
        int required = get_bit(ptr, 0) != 0;

        int value = xor_sum ^ required;

        assert(!get_bit(marks, base_col[eqn]));

        // should otpimize this
        propagated_lits.push_back(
            make_pair((value == 0 ? -1 : 1) * var_id[base_col[eqn] - 1], eqn));
        xup_level[eqn] = level;
        continue;
      }

      if (num_ones == 0)
        continue;

      unset_bit(remained_ones, base_col[eqn]);
      int next_watch = msb_position(remained_ones);
      set_watch(next_watch, eqn);
      assert(watch_col[eqn] != 0);
    }
  }

#ifndef NDEBUG
  sort(propagated_lits.begin(), propagated_lits.end());
  assert(unique(propagated_lits.begin(), propagated_lits.end()) ==
         propagated_lits.end());
#endif

  return propagated_lits;
}

inline void GaussianDense::rewatch(int col) {
  auto &w = watched_eqns[col];
  auto i = w.begin();
  auto end = w.end();

  for (auto j = i; j != end; ++j) {
    int eqn = *j;
    char *ptr = eqn_start_ptr(eqn);
    assert(watch_col[eqn] != col);
    if (watch_col[eqn] == 0 && base_col[eqn] != col && get_bit(ptr, col)) {
      watch_col[eqn] = col;
      *i++ = eqn;
    }
  }

  w.resize(i - w.begin());
  if (w.size() != 0)
    set_bit(watch_mask, col);
}

void GaussianDense::backtrack(int level) {
  auto i = assigned.rbegin();

  for (; i != assigned.rend(); ++i) {
    Var &v = internal->var(*i);
    if (v.level > level) {
      LOG("backtrack %d\n", *i);
      int pos = internal->get_gaussian_pos(*i);
      unset_bit(marks, pos);
      unset_bit(assignment, pos);
      rewatch(pos);
    } else {
      break;
    }
  }

  assigned.resize(assigned.rend() - i);
}

Clause *GaussianDense::conflict(int eqn) {
  START(xanalyze);
  char *ptr = eqn_start_ptr(eqn);
  int size = count(ptr);

  int bytes = Clause::bytes(size);
  Clause *reason = (Clause *)new char[bytes];
  reason->size = size;

  // FIXME: Still need a better implementation
  int j = 0;
  while (count(ptr)) {
    int msb = msb_position(ptr);
    assert(internal->val(var_id[msb - 1]) != 0);
    reason->literals[j++] = -internal->val(var_id[msb - 1]) * var_id[msb - 1];
    unset_bit(ptr, msb);
  }
  for (int i = 0; i < size; ++i)
    set_bit(ptr, internal->get_gaussian_pos(abs(reason->literals[i])));
  
  LOG("gen xconflict");
  
  STOP(xanalyze);
  return reason;
}

pair<Clause *, int> GaussianDense::analyze(int eqn) {
  START(xanalyze);
  LOG("Analyze conflict in block %d's %d equation", index, eqn);

  char *ptr = eqn_start_ptr(eqn);
  int size = count(ptr);

  int bytes = Clause::bytes(size);
  Clause *reason = (Clause *)new char[bytes];
  reason->size = size;

  // FIXME: Still needs a better implementation
  LOG("analyze BASE: %d", var_id[base_col[eqn] - 1]);

  assert(get_bit(ptr, base_col[eqn]));

  int j = 0;
  while (count(ptr)) {
    int msb = msb_position(ptr);
    assert(internal->val(var_id[msb - 1]) != 0);
    if (msb == base_col[eqn])
      reason->literals[j++] = internal->val(var_id[msb - 1]) * var_id[msb - 1];
    else
      reason->literals[j++] = -internal->val(var_id[msb - 1]) * var_id[msb - 1];
    unset_bit(ptr, msb);
  }
  for (int i = 0; i < size; ++i)
    set_bit(ptr, internal->get_gaussian_pos(abs(reason->literals[i])));

  STOP(xanalyze);
  return make_pair(reason, xup_level[eqn]);
}

bool GaussianDense::all_satisfied() {
  for (int i = 0; i < num_eqns; ++i) {
    char *ptr = eqn_start_ptr(i);
    int flag = 0;
    for (int j = 1; j <= num_vars; ++j) {
      if (get_bit(ptr, j) && internal->val(var_id[j - 1]) == 0) {
        assert(0 && "no all assigned");
        return false;
      }
      if (get_bit(ptr, j))
        flag ^= internal->val(var_id[j - 1]) > 0;
    }
    if (flag != !!get_bit(ptr, 0)) {
      LOG("ERROR: xclause not satisfied");
      assert(0 && "value is incorrect");
      return false;
    }
  }
  return true;
}

} // namespace CaDiCaL
