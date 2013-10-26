#include "prepare_states_bath.hpp"
#include <vector>
#include <map>
#include <tuple>
#include <iostream>
#include <fstream>
#include <random>
#include <mpi.h>
#include "make_wb.hpp"
#include "make_cb.hpp"
#include "make_bath.hpp"
#include "make_initial_state.hpp"
//#include "gnuplot_i.hpp"
#include "o.hpp"

std::tuple<
  std::vector<double> ,
  std::vector<double> ,
  std::vector<std::vector<double>> 
>
prepare_states_bath(
  std::map<std::string,double> sp,
  std::map<std::string,double> mp,
  std::map<std::string,double> ip
    ) {
  using namespace std;
  //cout << "Dealing with bath ... \n";
  vector<double> wb,cb;
  wb = make_wb(sp["n_modes"], sp["dw"]);  o(wb,"wb.dat"); 
  cb = make_cb(wb,sp["eta"],sp["w_c"],sp["dw"]); o(cb,"cb.dat"); 

  //cout << "Making initial states ... \n";
  vector<vector<double>> states(mp.at("n_trajs"));
  int  seed = mp.at("seed");
  int world_size; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  //std::random_device rd;
  //std::mt19937 gen(rd());
  //std::mt19937 gen(world_rank);
  std::mt19937 gen(world_rank+seed);
  std::uniform_int_distribution<> dis(0,std::numeric_limits<int>::max());
  seed=dis(gen);
  for(auto &i:states) {
    vector<double> pb,qb; tie(pb,qb) = make_bath(wb,cb,sp.at("beta"),++seed);
    i= make_initial_state(ip.at("pp"),ip.at("qp"),ip.at("n1"), mp.at("bin_start"),mp.at("n_shift"),pb,qb,++seed); 
  }
  string ost = "states/"+to_string(world_rank)+".dat";
  o(states,ost);
  //o(states);

  //Gnuplot g;
  //vector<double> x,y;
  //for(auto &i:states) {
    //x.push_back(i.at(7+wb.size()+0));
    //y.push_back(i.at(7+wb.size()+1));
  //}
  //g.cmd("set terminal x11 persist");
  //g.plot_xy(x,y);
  //g.cmd("pause 1");
  //std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );

  return make_tuple(wb,cb,states);
}
