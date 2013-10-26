#include "traj_step.hpp"

#include <vector>
#include <map>
#include <functional>
//#include <tuple>
#include <algorithm>
#include <iostream>
#include <boost/numeric/odeint.hpp>

void traj_step::operator()(std::vector<double> &state) {
  using namespace boost::numeric::odeint;
  bulirsch_stoer< std::vector<double> > stepper(1e-6,0.0,0.0,0.0);
  typedef std::vector<double> state_type ;
  adams_bashforth_moulton< 5 , state_type > abm;

  //integrate_const(abm,eom, state, state[0],state[0]+dt,dt*1e-4);
  integrate(eom, state, state[0],state[0]+dt,dt*1e-3);
  state[0]+=dt;
  state.back()=ham(state);
}
