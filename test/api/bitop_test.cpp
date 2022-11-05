#include "../../src/bitops.hpp"
#include <iostream>
#include <vector>
#include <cassert>

#define size 8
using namespace std;

int main() {
  char a[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  char b[8] = {7, 6, 5, 4, 3, 2, 1, 0};
  char c[8] = {1, 3, 5, 7, 9, 2, 4, 6};
  char d[8];
  
  Bits A(size), B(size);
  A.load_from_ptr(a);
  B.load_from_ptr(b);
  {
    Bits C = A & B;
    
    for (int i = 0; i < 8; ++i)
      assert(C.arr[i] == (a[i] & b[i]));
    C.store_to_ptr(d);
    for (int i = 0; i < 8; ++i)
      assert(d[i] == (a[i] & b[i]));
  }
}
