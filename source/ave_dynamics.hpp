#pragma once

#include <vector>
#include <map>
#include <tuple>

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
    ) ;
