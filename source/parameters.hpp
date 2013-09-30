#pragma once
#include <tuple>
#include <map>

std::tuple<
  std::map<std::string,double>,
  std::map<std::string,double>,
  std::map<std::string,double>,
  std::map<std::string,double>
> 
get_param(const char *lua_config_filename) ;
