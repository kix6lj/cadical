#ifndef _gaussian_hpp_INCLUDE
#define _gaussian_hpp_INCLUDE

#include <bitset>

#include "inttypes.hpp"
#include <vector>

namespace CaDiCaL {
struct Internal;

struct XOR_Clause {};

// This module is reponsible for handling a connected block of XOR equations
//
class Gaussian {
  vector<XOR_Clause *> equations;

public:
  Gaussian();

  // Deduce a number of assignments or conflict from some newly assigne
  // variables
  //
  int deduce(vector<int> assigned_vars);

  // When performing conflict analysis, analyze the root cause of some
  // deduced variables
  //
  int analyze_cause(vector<int> deduced_vars);
};

struct Raw_XOR_Clause {
  int flip;
  vector<int> vars;
};

struct Raw_XOR_Equations {
  vector<Raw_XOR_Clause *> equations; // raw xor-equations
  vector<int> parent;         // parent of equation in the union-find set
  vector<int> component_size; // size of each union-find set
  Internal *internal;
  
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
  void analyze();
};
}; // namespace CaDiCaL

#endif
