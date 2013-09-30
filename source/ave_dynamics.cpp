#include "ave_dynamics.hpp"

#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

template <typename... T>
auto zip_range(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
    auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}

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
  vvd ave_traj, ave_n;
  vector<double>  time;
  double dt = tp["time_step"];
  copy(irange(0,static_cast<int>(tp["n_times"])) | transformed( [dt](double itime ){return itime*dt;}),back_inserter(time));
  int nt = tp["n_times"];
  vector<vector<double>> ave_observes,ave_state(time.size()),ave_nn(time.size());
  transform(
      time, 
      make_zip_iterator(
        boost::make_tuple(begin(ave_state),begin(ave_nn))
        ), 
      [](double t) {return boost::make_tuple(vd(2,t),vd(2,0));});
  for(auto i:time) cout << i << "   " ;
  cout <<"\n";
  for(auto i:ave_observes) { for(auto j:i) { cout << j<< "  ";} cout << "\n";}
  cout <<"\n";
  return std::make_tuple(ave_traj,ave_n);
}
