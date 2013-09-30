#include "parameters.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lua_read.hpp"

std::tuple<
  std::map<std::string,double>,
  std::map<std::string,double>,
  std::map<std::string,double>,
  std::map<std::string,double>
> 
get_param(const char *lua_config_filename) {
using namespace std;
  std::vector<std::map<std::string,double>> ret;

  ofstream ofpara("param.dat");

  auto sp = Lua_read(lua_config_filename).read_config(vector<string>{"m","e1","e2","w1","w2","rp1","rp2","v12","eta","hbar","w_c","dw","n_modes","n_shift","boltzmann_constant","beta"});
  ofpara << "map of system's parameters: \n"; for(auto i: sp) ofpara<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto ip = Lua_read(lua_config_filename).read_config(vector<string>{"pp","qp","n1","n2"});
  ofpara << "map of initial parameters: \n"; for(auto i: ip) ofpara<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto tp = Lua_read(lua_config_filename).read_config(vector<string>{"n_times","end_time","time_step"});
  ofpara << "map of integration parameters: \n"; for(auto i: tp) ofpara<< setw(20) << i.first << " = " << i.second << "\n";
  
  auto mp = Lua_read(lua_config_filename).read_config(vector<string>{"seed","n_shift","bin_start","bin_end","n_trajs"});
  ofpara << "map of method parameters: \n"; for(auto i: mp) ofpara<< setw(20) << i.first << " = " << i.second << "\n";
  
  return make_tuple(sp,ip,tp,mp);
}
