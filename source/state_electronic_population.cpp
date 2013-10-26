#include "state_electronic_population.hpp"
#include <vector>
#include <map>
//#include <tuple>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>

#include <mpi.h>

#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
//#include <boost/assign/std/vector.hpp>
//#include <boost/phoenix/phoenix.hpp>

#include "o.hpp"
#include "state.hpp"
#include "box.hpp"

std::vector<double> state_electronic_population( std::vector<std::vector<double>> states, double n1, double n2, double bin_end, double n_shift) {
  //auto el_pop = electronic_population(trajs,ip.at("n1"),ip.at("n2"),mp.at("bin_end"),mp.at("n_shift"));
  using std::cout;
  using std::vector;
  using namespace boost;
  using namespace boost::adaptors;
  typedef vector<double> vd; typedef vector<vector<double>> vvd; typedef vector<vector<vector<double>>> vvvd;
  //vvd nnps;
  //boost::transform(states,back_inserter(nnps),[] (vd state) {
      //return state;
      //});
  vd nnp = std::accumulate(states.begin(),states.end(),vd(4,0.0),[=](vd nnp_prev,vd s)->vd {
      auto time = s.front();
      auto n1b = box(get_n{n_shift}(s,1), n1, bin_end);
      auto n2b = box(get_n{n_shift}(s,2), n1, bin_end);
      auto nnb = n1b-n2b;
      vd cnnp ;
      cnnp.push_back(time); 
      cnnp.push_back(n1b); 
      cnnp.push_back(n2b); 
      cnnp.push_back(nnb); 
      //vd cnnp {time, n1b, n2b, nnb};
      std::transform(nnp_prev.begin()+1,nnp_prev.end(),cnnp.begin()+1,cnnp.begin()+1,[](double a, double b) { return b+a;});
      //boost::for_each(nnp_prev,cnnp,[](double a, double &b) { b+=a;});
      return cnnp;
      });
  //for(auto &i:nnp) i/=states.size();
  std::for_each(nnp.begin()+1,nnp.end(),[&states](double &i){i/=states.size();});
  //using namespace boost::phoenix;
  //using namespace boost::phoenix::placeholders;

  //auto sts = states.size();
  //for_each(nnp.begin()+1,nnp.end(),lambda[arg1/=states.size()])(arg1,sts);
  //for_each(nnp.begin()+1,nnp.end(),lambda(ss=arg2)[arg1/=states.size()])(arg1,sts);
  int world_size; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  string nnst  = "ele_nn/"+to_string(world_rank)+".dat";
  o(nnp,nnst);
  string ost = "ele_states/"+to_string(world_rank)+".dat";
  o(states,ost);

  vector<double> sub_avgs;
  //if (world_rank == 0) {
    sub_avgs.resize(world_size*nnp.size());
  //}
  //MPI_Allgather(&nnp[0], 1, MPI_DOUBLE, &sub_avgs[0], 1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(nnp.data(), nnp.size(), MPI_DOUBLE, sub_avgs.data(), nnp.size(), MPI_DOUBLE, MPI_COMM_WORLD);
  //MPI_Gather(nnp.data(), 1, MPI_DOUBLE, sub_avgs.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //if(world_rank==1) {
    //cout << "nnp 1  ";
    //for(auto i: nnp) cout << i << "  ";
    //cout << "\n";
  //}
  //if(world_rank==0) {
    //cout << "nnp 0  ";
    //for(auto i: nnp) cout << i << "  ";
    //cout << "\n";
    //cout << "sub_avgs ";
    //for(auto i: sub_avgs) cout << i << "  ";
    //cout << "\n";
  //}

  //if (world_rank == 0) {
    vector<vector<double>> mpi_nnps;
    //vector<vector<double>> mpi_nnps(world_size);
    for (size_t i = 0; i < world_size; ++i) {
      //mpi_nnps.push_back(vector<double>(nnp.size(),0.1));
      mpi_nnps.push_back(vector<double>(sub_avgs.begin()+i*nnp.size(),sub_avgs.begin()+(i+1)*nnp.size()));
    }
    //vd mpi_nnp = mpi_nnps.front();
    vd mpi_nnp = std::accumulate(mpi_nnps.begin(), mpi_nnps.end(),vd(4,0.0),[=](vd a, vd b)-> vd {
        vd c(a);
        std::transform(a.begin()+1,a.end(),b.begin()+1, c.begin()+1, [] (double a, double b) {return a+b;});
        return c;
        });

  //if(world_rank==0) {cout << "-- nnp "; std::copy(nnp.begin(),nnp.end(),ostream_iterator<double>(cout,"  ")); cout << "\n";}
  //if(world_rank==0) {cout << "-- mpi_nnp "; std::copy(mpi_nnp.begin(),mpi_nnp.end(),ostream_iterator<double>(cout,"  ")); cout << "\n";}
  return mpi_nnp;
  //return vector<double>(nnp.size(),0.1);
  //return nnp;
}

