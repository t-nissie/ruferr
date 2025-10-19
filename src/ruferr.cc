// ruferr.cc
////
#include <stdexcept> // for std::invalid_argument
#include <iostream>
#include "coordinates_class.h"

int main(int argc, char *argv[])
{
  int howmany, Lx, Ly, Lz, padding_y, n_times;
  if (argc==7) {
    howmany   = atoi(argv[1]);
    Lx        = atoi(argv[2]);
    Ly        = atoi(argv[3]);
    Lz        = atoi(argv[4]);
    padding_y = atoi(argv[5]);
    n_times   = atoi(argv[6]);
  } else {
    throw std::invalid_argument("The number of command line arguments should be 6.");
  }
  Coordinates coords(howmany, Lx, Ly, Lz, padding_y);
  coords.set_values();
  coords.forward_backward(n_times);
  coords.show_some();

  return 0;
}
// local variables:
//   compile-command: "make clean; make -j4 && ./ruferr 3 32 64 128 3 100"
// End:
