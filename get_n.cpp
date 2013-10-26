#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>

using namespace std;

int main(int argc, const char *argv[])
{
  vector<vector<double>> px;
  ifstream inf(argv[1]);

  for(string line;getline(inf, line);) if(!inf.eof()) {
    istringstream ll(line);
    px.emplace_back((istream_iterator<double>(ll)), istream_iterator<double>());
  }

  //for(string line;getline(inf, line);) {
    //istringstream iss(line);
    //vector<double> l;
    //double x;
    //if(!inf.eof()) {
      //while(iss>> x) l.push_back(x);
      //px.push_back(l);
    //}
  //}

  ofstream of("ns.dat");
  cout << "size = " << px.size() << "\n";
  for(auto i:px) { 
    vector<double> nn = {i[0], 0.5*(i[3]*i[3]+i[4]*i[4])-0.35, 0.5*(i[5]*i[5]+i[6]*i[6])-0.35};
    boost::copy(nn, ostream_iterator<double>(of,"  ")); 
    of << "\n";
  }
  return 0;
}
