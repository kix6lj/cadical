#ifndef _visualization_hpp_INCLUDE
#define _visualization_hpp_INCLUDE

namespace CaDiCaL {
struct Internal;

struct GraphVisualizer {
  Internal *internal;
  GraphVisualizer(Internal *internal) : internal(internal) {}
  void build_from_cnf();
  void write_graph_dot(FILE *file);
};

} // namespace CaDiCaL
#endif
