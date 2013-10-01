#include "ave_dynamics.hpp"

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

#include <boost/numeric/odeint.hpp>

#include "pcet.hpp"
#include "traj_step.hpp"
#include "ave_over_states.hpp"


std::vector<double> el_pop( std::vector<std::vector<double>> states) {
  return states.front();
}

void trajs_step( std::vector<std::vector<double>> &states, 
    function<void(std::vector<double>&)> my_traj_step
    ) {
  boost::for_each(states,my_traj_step);
}

struct avs {
  std::vector<std::vector<double>> states;
  function<void(std::vector<double>&)> my_traj_step;
  boost::tuple<std::vector<double>,std::vector<double>> operator()(double t) {
    std::vector<double> ave_state = ave_over_states(states);
    std::vector<double> ave_nn = el_pop(states);
    trajs_step(states, my_traj_step); ///> step
    return boost::make_tuple(ave_state,ave_nn);
  }
};

std::tuple<
  std::vector<std::vector<double>> ,
  std::vector<std::vector<double>> 
>
ave_dynamics(
  std::vector<double> wb,
  std::vector<double> cb,
  std::vector<std::vector<double>> states,
  std::map<std::string,double> sp,
  std::map<std::string,double> tp
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

  vvd ave_traj, ave_n;
  vector<double>  time;
  copy(irange(0,static_cast<int>(tp["n_times"])) | transformed( [dt](double itime ){return itime*dt;}),back_inserter(time));
  int nt = tp["n_times"];
  vector<vector<double>> ave_observes,ave_state(time.size()),ave_nn(time.size());
  transform( time, make_zip_iterator( boost::make_tuple(begin(ave_state),begin(ave_nn))), 
      avs{states,my_traj_step}
      );
      //[](double t) {
      //return boost::make_tuple(vd(2,t),vd(2,0));});
  for(auto i:time) cout << i << "   " ;
  cout <<"\n";
  for(auto i:ave_observes) { for(auto j:i) { cout << j<< "  ";} cout << "\n";}
  cout <<"\n";
  return std::make_tuple(ave_traj,ave_n);
}
