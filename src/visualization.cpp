#include "internal.hpp"
#include "visualization.hpp"

namespace CaDiCaL {
  void GraphVisualizer::build_from_cnf() {

  }
  void GraphVisualizer::write_graph_dot(FILE *file) {
    fprintf(file, "graph cnf_graph {\n");
    // variables
    fprintf(file, "  node [shape=point];\n");
    for (auto idx : internal->vars)
      fprintf(file, "  l%d [color=red];\n", idx);
    for (int i = 0; i < internal->clauses.size(); ++i) {
      fprintf(file, "  c%d [color=blue];\n", i);
      auto c = internal->clauses[i];
      for (int j = 0; j < c->size; ++j)
	fprintf(file, "  c%d -- l%d\n", i, abs(c->literals[j]));
    }
    fprintf(file, "}\n");
  }
};
