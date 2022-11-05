#ifndef _gaussian_hpp_INCLUDE
#define _gaussian_hpp_INCLUDE

#include <bitset>

#include "inttypes.hpp"
#include <vector>

namespace CaDiCaL {
struct Internal;
struct Clause;
struct Raw_XOR_Equations;

struct Raw_XOR_Clause {
  int flip;
  vector<int> vars;
};

class GaussianDense{
private:
  Internal *internal;
  Raw_XOR_Equations *manager;
  int index;
  int num_bytes_per_eqn;
  int num_eqns;
  int num_vars;
  char *equations;  // all clauses
  char *assignment; // assignment of variables
  char *marks;      // is a variable assigned?
  char *base_mask;
  char *watch_mask;
  char *temp_mem;
  vector<int> var_id;
  vector<int> base_col;             // base-column of each equation
  vector<int> base_row;             // row of each base column
  vector<int> watch_col;            // watch-column of each equation
  vector<char> col_is_base;         // Is the column a base?
  vector<char> col_is_watch;        // Is the column a watch?
  vector<vector<int>> watched_eqns; // Watched rows of each column
  vector<int> assigned;             // assigned variables from the solver
  vector<int> xup_level;
  
  /* =----------------------------------= */
  // Internal Functions
  //
  inline char *eqn_start_ptr(int i) {
    return equations + num_bytes_per_eqn * i;
  }

  // *(i + k) ^= *(j + k)
  // *(r + k) = *(i + k) ^ *(j + k)
  void XOR(char *i, char *j);

  void XOR(char *r, char *i, char *j);

  void OR(char *i, char *j);

  void OR(char *r, char *i, char *j);

  // *(r + k) = *(i + k) & ~*(j + k)
  void IAND(char *i, char *j);

  void IAND(char *r, char *i, char *j);

  void AND(char *i, char *j);

  void AND(char *r, char *i, char *j);

  int msb_position(char *ptr);

  int lsb_position(char *ptr);

  char get_bit(char *ptr, int i);

  void set_bit(char *ptr, int i);

  void unset_bit(char *ptr, int i);

  int count(char *ptr);

  void set_watch(int col, int eqn);

  void remove_watch(int col, int eqn);

  void rewatch(int col);

  void init_watch();

  void print(char *ptr);
  
  Clause * conflict(int eqn);
public:
  GaussianDense(Internal *internal, Raw_XOR_Equations *manager, int index);
  ~GaussianDense();

  void initialize(const vector<Raw_XOR_Clause *> &original);

  vector<pair<int, int>> propagate(const vector<int> &var,
                                   const vector<int> &assgns,
				   int level);
  void backtrack(int level);

  pair<Clause*, int> analyze(int eqn);

  bool all_satisfied();
};

struct Raw_XOR_Equations {
  vector<Raw_XOR_Clause *> equations; // raw xor-equations
  vector<int> parent;         // parent of equation in the union-find set
  vector<int> component_size; // size of each union-find set
  vector<GaussianDense *> blocks;
  Internal *internal;
  bool conflict_flag;
  Clause *conflict;
  vector<pair<Clause*, int>> generated_reasons;
  
  Raw_XOR_Equations(Internal *internal);
  ~Raw_XOR_Equations();

  void add_clause(const vector<int> &lits, int flip);

  // utilities of unio-find set
  //
  void merge(int x, int y);
  int find(int x);

  // check if there are same variable in two clauses
  // using the mark array in the internal solver
  //
  bool check_incidence(const vector<int> &vars1, const vector<int> &vars2);
  void analyze_connectivity();
  void partition();
  void analyze();

  /* =--------------------------------------= */
  // Functions called by the main loop
  /* =--------------------------------------= */

  // distribute the propagated vars to each block
  vector<pair<int, int>> propagate(const vector<int> &lit, int level);

  // call backtrack to all blocks
  void backtrack(int level);

  Clause *lazily_gen_reason(int block, int eqn);

  bool all_satisfied();
};
}; // namespace CaDiCaL

#endif
