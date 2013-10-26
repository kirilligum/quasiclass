#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <limits>
#include <tuple>
#include <mpi.h>
//#include <boost/range/algorithm.hpp>
//#include <boost/range/algorithm_ext.hpp>
//#include <boost/range/numeric.hpp>
//#include <boost/math/constants/constants.hpp>
//#include <boost/numeric/odeint.hpp>

#include "o.hpp"
#include "state.hpp"
#include "integrate_traj.hpp"
#include "box.hpp"
#include "electronic_population.hpp"
#include "parameters.hpp"
//#include "gnuplot_i.hpp"
#include "prepare_states_bath.hpp"
#include "ave_dynamics.hpp"

using namespace std;
using namespace std::placeholders;

int main(int argc, const char *argv[]) {
  MPI_Init(nullptr,nullptr);
  using namespace std; 
  typedef vector<double> vd; typedef vector<vector<double>> vvd; typedef vector<vector<vector<double>>> vvvd;

  //cout << "Reading parameters from " << argv[1] << " into the pro and into param.dat ... \n";
  map<string,double> sp,ip,tp,mp;
  tie(sp,ip,tp,mp)=get_param(argv[1]);

  //cout << "Preparing initial states and bath ... \n";
  vector<double> wb,cb;
  vector<vector<double>> states(mp.at("n_trajs"));
  tie(wb,cb,states) = prepare_states_bath(sp,mp,ip);

  vvd ave_traj, ave_n;
  tie(ave_traj,ave_n) = ave_dynamics(wb,cb,states,sp,tp,ip,mp);

  //cout << "Integrating trajectories \n";
  //vector<vector<vector<double>>> trajs;
  //transform(states.begin(),states.end(),back_inserter(trajs), bind(integrate_traj,_1,tp["n_times"],tp["end_time"],sp,wb,cb));

  //auto traj = trajs.front();
  //for(auto &i:traj) {
    //i.push_back(get_n{mp.at("n_shift")}(i,1));
    //i.push_back(get_n{mp.at("n_shift")}(i,2));
  //}
  //o(traj,"traj.dat");

  //cout << "binning ... \n";
  //auto el_pop = electronic_population(trajs,ip.at("n1"),ip.at("n2"),mp.at("bin_end"),mp.at("n_shift"));
  //o(el_pop ,"el_pop.dat");
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
