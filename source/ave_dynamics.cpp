#include "ave_dynamics.hpp"

#include <vector>
#include <map>
//#include <tuple>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

#include <boost/numeric/odeint.hpp>

#include "pcet.hpp"
#include "traj_step.hpp"
#include "ave_dyn_step.hpp"

#include "o.hpp"

std::tuple<
  std::vector<std::vector<double>> ,
  std::vector<std::vector<double>> 
>
ave_dynamics(
  std::vector<double> wb,
  std::vector<double> cb,
  std::vector<std::vector<double>> states,
  std::map<std::string,double> sp,
  std::map<std::string,double> tp,
  std::map<std::string,double> ip,
  std::map<std::string,double> mp
    ) {
  //using namespace std;
  using std::cout;
  using std::begin;
  using std::end;
  using std::vector;
  using namespace boost;
  using namespace boost::adaptors;
  typedef vector<double> vd; typedef vector<vector<double>> vvd; typedef vector<vector<vector<double>>> vvvd;

  auto my_energy = pcet(sp,wb,cb);
  auto my_eom = eom(d_pcet(sp,wb,cb));
  double dt = tp["time_step"];
  auto my_traj_step = traj_step{dt,my_energy,my_eom};

  rmkdir("int_states");
  rmkdir("energies");
  vector<double>  time;
  copy(irange(0,static_cast<int>(tp["n_times"])) | transformed( [dt](double itime ){return itime*dt;}),back_inserter(time)); ///> make a vector of times
  vvd ave_traj(time.size()), ave_n(time.size());
  transform( time, make_zip_iterator( boost::make_tuple(begin(ave_traj),begin(ave_n))), 
      ave_dyn_step{states,my_energy,my_traj_step,ip,mp}
      );

  ofstream ofan("ave_n.dat");
  for(auto i:ave_n) {
    for(auto j:i) ofan << j << "   " ;
    ofan <<"\n";
  }

  ofstream ofat("ave_traj.dat");
  for(auto i:ave_traj) {
    for(auto j:i) ofat << j << "   " ;
    ofat <<"\n";
  }

  ofstream oftime("time.dat");
  for(auto i:time) oftime << i << "   " ;
  oftime <<"\n";
  
  return std::make_tuple(ave_traj,ave_n);
}
