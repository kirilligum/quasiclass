#pragma once

#include <vector>
#include <map>
#include <boost/tuple/tuple.hpp>

struct ave_dyn_step {
  std::vector<std::vector<double>> states;
  std::function<double(std::vector<double>&)> my_energy;
  std::function<void(std::vector<double>&)> my_traj_step;
  std::map<std::string,double> ip;
  std::map<std::string,double> mp;
  boost::tuple<std::vector<double>,std::vector<double>> operator()(double t);
};

