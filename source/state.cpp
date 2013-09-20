#include "state.hpp"

#include <vector>
#include <iostream>
using namespace std;

double get_n::operator()(vector<double> s,int i) {
  if(i==1 || i==2) {
    auto j=i;
    j*=2;
    j+=1;
    return 0.5*(s[j]*s[j]+s[j+1]*s[j+1])-n_shift;
  } else {
    cout << " error wrong N in nN \n";
    return 1337.35505;
  }
}

