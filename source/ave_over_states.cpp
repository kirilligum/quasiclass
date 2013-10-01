#include "ave_over_states.hpp"
#include <vector>
#include <map>
//#include <tuple>
#include <algorithm>
#include <iostream>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>



std::vector<double> ave_over_states( std::vector<std::vector<double>> states) {
  using std::cout;
  //std::cout << "averaging ... \n";
  //cout << "states.size() = " << states.size() << "\n";
  //cout << "states.front().size() = " << states.front().size() << "\n";
  //for(auto i:states) {
    //cout << "s  ";
    //for (auto j:i) cout << "  " << j;
    //cout << "\n";
  //}
  //std::vector<double> ret;
  //std::vector<double> ret(states.front().size(),0.0); ///> ask Steve.  why is this taking 40 secnds
  //for (size_t i = 0; i < states.size(); ++i) {
    //for (size_t j = 0; j < states.front().size(); ++j) {
      ////ret[j]+=states[i][j];
    //}
  //}
  //return std::accumulate(
  auto ret =  std::accumulate(
      states.begin()+1,
      states.end(),
      states.front(),
      [](std::vector<double> a,std::vector<double> b){ 
        //std::cout << " i   " ;
        std::vector<double> c; 
        //std::vector<double> c(a.size()); 
        boost::transform(a,b,back_inserter(c), [](double aa, double bb) {
          return aa+bb;}
        ); 
        return c;
      });
  for(auto i:ret) {i/=states.size();}
  //for(auto i:ret)  cout << i << "  "; cout << " <--- ave_state\n";
  return ret;
  //return states.front();
}
