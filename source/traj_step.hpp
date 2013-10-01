#pragma once

#include <vector>
#include <functional>

struct traj_step {
  //std::vector<double> state;
  const double dt;
  std::function<double(const std::vector<double>&)> ham;
  std::function<void(const std::vector<double>&,std::vector<double>&,const double)> eom;
  void operator()(std::vector<double> &state) ;
};

