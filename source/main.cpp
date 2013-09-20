#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <numeric>
#include <iomanip>
//#include <boost/range/algorithm.hpp>
//#include <boost/range/algorithm_ext.hpp>
//#include <boost/range/numeric.hpp>
//#include <boost/math/constants/constants.hpp>
//#include <boost/numeric/odeint.hpp>

#include "o.hpp"
#include "make_wb.hpp"
#include "make_cb.hpp"
#include "make_bath.hpp"
#include "make_initial_state.hpp"
#include "state.hpp"
#include "integrate_traj.hpp"
#include "box.hpp"
#include "electronic_population.hpp"
#include "lua_read.hpp"

using namespace std;
using namespace std::placeholders;

int main(int argc, const char *argv[]) {
  using namespace std; 
  typedef vector<double> vd; typedef vector<vector<double>> vvd; typedef vector<vector<vector<double>>> vvvd;
  cout << "Reading parameters from " << argv[1] << " ... \n";

  auto sp = Lua_read(argv[1]).read_config(vector<string>{"m","e1","e2","w1","w2","rp1","rp2","v12","eta","hbar","w_c","dw","n_modes","n_shift","boltzmann_constant","beta"});
  cout << "map of system's parameters: \n"; for(auto i: sp) cout<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto ip = Lua_read(argv[1]).read_config(vector<string>{"pp","qp","n1","n2"});
  cout << "map of initial parameters: \n"; for(auto i: ip) cout<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto tp = Lua_read(argv[1]).read_config(vector<string>{"n_times","end_time","time_step"});
  cout << "map of integration parameters: \n"; for(auto i: tp) cout<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto mp = Lua_read(argv[1]).read_config(vector<string>{"seed","n_shift","bin_start","bin_end","n_trajs"});
  cout << "map of method parameters: \n"; for(auto i: mp) cout<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto seed = mp.at("seed");


  cout << "Dealing with bath ... \n";
  auto wb = make_wb(sp["n_modes"], sp["dw"]);  o(wb,"wb.dat"); 
  auto cb = make_cb(wb,sp["eta"],sp["w_c"],sp["dw"]); o(cb,"cb.dat"); 

  cout << "Making initial states ... \n";
  vector<vector<double>> states(mp.at("n_trajs"));
  for(auto &i:states) {
    vector<double> pb,qb; tie(pb,qb) = make_bath(wb,cb,sp.at("beta"),++seed);
    i= make_initial_state(ip.at("pp"),ip.at("qp"),ip.at("n1"), mp.at("bin_start"),mp.at("n_shift"),pb,qb,++seed); 
  }
  o(states,"states.dat");

  
  //vector<double> dich;
  //transform(states.begin(),states.end(),back_inserter(dich),[](vector<double>v){return (v[1]);});
  //sort(dich.begin(),dich.end());
  //ofstream of("dich.dat");
  //for(auto i:dich) of<<i<<"\n";

  cout << "Integrating trajectories \n";
  vector<vector<vector<double>>> trajs;
  transform(states.begin(),states.end(),back_inserter(trajs), bind(integrate_traj,_1,tp["n_times"],tp["end_time"],sp,wb,cb));

  double an1=0, an2=0;
  for(auto is: trajs) {
    an1+=is.front().at(207);
    an2+=is.front().at(208);
  }
  cout << "tot n = " << an1 << "  " << an2 << "\n";
  cout << "ave n = " << an1/states.size() << "  " << an2/states.size() << "\n";
  
  auto traj = trajs.front();
  for(auto &i:traj) {
    i.push_back(get_n{mp.at("n_shift")}(i,1));
    i.push_back(get_n{mp.at("n_shift")}(i,2));
  }
  o(traj,"traj.dat");

  cout << "binning ... \n";
  auto el_pop = electronic_population(trajs,ip.at("n1"),ip.at("n2"),mp.at("bin_end"),mp.at("n_shift"));
  o(el_pop ,"el_pop.dat");
  
  return 0;
}
