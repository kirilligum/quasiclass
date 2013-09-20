#include "tunnel.hpp"

vector<double> tunnel::make_initial_state(const map<string,double> param){
  vector<double> ret;
  return ret;
}

double tunnel::energy(vector<double> state) {
  return kin+v1*n1+v2*n2+v12*nn;
  return 0.0;
}

void tunnel::eom(const vector<double> &state, vector<double> &dstates, const double) {
}
