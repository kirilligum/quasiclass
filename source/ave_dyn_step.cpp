#include "ave_dyn_step.hpp"

#include <vector>
#include <map>
//#include <tuple>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

//#include <boost/numeric/odeint.hpp>

#include "pcet.hpp"
#include "traj_step.hpp"
#include "ave_over_states.hpp"
#include "state_electronic_population.hpp"

#include "o.hpp"

void trajs_step( std::vector<std::vector<double>> &states, 
    function<void(std::vector<double>&)> my_traj_step
    ) {
  boost::for_each(states,my_traj_step);
}

boost::tuple<std::vector<double>,std::vector<double>> ave_dyn_step::operator()(double t) {
  int world_size; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  string ost = "int_states/"+to_string(world_rank)+".dat";
  o(states,ost);
  std::vector<double> ave_state = ave_over_states(states);
  std::vector<double> ave_nn = state_electronic_population(states,ip.at("n1"),ip.at("n2"),mp.at("bin_end"),mp.at("n_shift"));
  trajs_step(states, my_traj_step); ///> step
  return boost::make_tuple(ave_state,ave_nn);
}

