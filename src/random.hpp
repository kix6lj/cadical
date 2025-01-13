#ifndef _random_hpp_INCLUDED
#define _random_hpp_INCLUDED

#include <cstdint>

// Random number generator.

namespace CaDiCaL {

/* Minimal PCG Generator */
struct pcg32_random_t { 
  uint64_t state;  
  uint64_t inc;
  
  // FIXME: dumb function to pass compilation
  operator uint64_t() const { return state; } 
 };

uint32_t pcg32_random_r(pcg32_random_t* rng);

class Random {

  pcg32_random_t rng;

  void add (uint64_t a) {
    if (!(rng.state += a))
      rng.state = 1;
    rng.inc += a;
    pcg32_random_r(&rng);
  }

public:
  // Without argument use a machine, process and time dependent seed.
  //
  Random ();

  Random (uint64_t seed) : rng ({seed, seed<<1}) { pcg32_random_r(&rng); }
  void operator= (pcg32_random_t seed) { rng = seed; }
  Random (const Random &other) : rng (other.seed ()) {}

  void operator+= (uint64_t a) { add (a); }
  pcg32_random_t seed () const { return rng; }

  uint64_t next() {
    // FIXME: Temporary solution to pass the compilation
    return (uint64_t) generate() * generate();
  }

  uint32_t generate () {
    return pcg32_random_r(&rng);
  }

  int generate_int () { return (int) generate (); }
  bool generate_bool () { return generate () < 2147483648u; }

  // Generate 'double' value in the range '[0,1]' excluding '1'.
  //
  double generate_double () { return generate () / 4294967295.0; }

  // Generate 'int' value in the range '[l,r]'.
  //
  int pick_int (int l, int r) {
    assert (l <= r);
    const unsigned delta = 1 + r - (unsigned) l;
    unsigned tmp = generate (), scaled;
    if (delta) {
      const double fraction = tmp / 4294967296.0;
      scaled = delta * fraction;
    } else
      scaled = tmp;
    const int res = scaled + l;
    assert (l <= res);
    assert (res <= r);
    return res;
  }

  int pick_log (int l, int r) {
    assert (l <= r);
    const unsigned delta = 1 + r - (unsigned) l;
    int log_delta = delta ? 0 : 32;
    while (log_delta < 32 && (1u << log_delta) < delta)
      log_delta++;
    const int log_res = pick_int (0, log_delta);
    unsigned tmp = generate ();
    if (log_res < 32)
      tmp &= (1u << log_res) - 1;
    if (delta)
      tmp %= delta;
    const int res = l + tmp;
    assert (l <= res), assert (res <= r);
    return res;
  }

  // Generate 'double' value in the range '[l,r]'.
  //
  double pick_double (double l, double r) {
    assert (l <= r);
    double res = (r - l) * generate_double ();
    res += l;
    assert (l <= res);
    assert (res <= r);
    return res;
  }
};

} // namespace CaDiCaL

#endif
