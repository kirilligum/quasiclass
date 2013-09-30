#pragma once
#include <vector>
#include <map>
#include <tuple>

std::tuple<
  std::vector<double> ,
  std::vector<double> ,
  std::vector<std::vector<double>> 
>
prepare_states_bath(
  std::map<std::string,double> sp,
  std::map<std::string,double> mp,
  std::map<std::string,double> ip
    ) ;
