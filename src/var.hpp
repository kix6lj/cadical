#ifndef _var_hpp_INCLUDED
#define _var_hpp_INCLUDED

namespace CaDiCaL {

struct Clause;

// This structure captures data associated with an assigned variable.

struct XOR_Reason {
  int block;
  int eqn;
};
  
struct Var {

  // Note that none of these members is valid unless the variable is
  // assigned.  During unassigning a variable we do not reset it.

  int level; // decision level
  int trail; // trail height at assignment
  bool from_xor;
  bool reason_flag;
  union {
    Clause *reason; // implication graph edge during search
    XOR_Reason xor_reason;
  };
};

} // namespace CaDiCaL

#endif
