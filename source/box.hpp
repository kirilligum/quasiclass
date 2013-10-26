#pragma once
#include <vector>

inline double box(double nt, double n, double half_bin) {
  if(abs(nt-n)<=half_bin) return 1.0;
  //if(abs(nt-n)<=half_bin) return 0.5/half_bin;
  else return 0.0;
}
