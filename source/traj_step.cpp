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
    integrate(eom, state, state[0],state[0]+dt,dt*1e-3);
    state[0]+=dt;
    state.back()=ham(state);
}
