#include "o.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/range/algorithm.hpp>
#include <boost/filesystem.hpp>

using namespace std;

void rmkdir( string dirname) {
  using namespace boost::filesystem;
  remove_all(dirname+"/*");
  if (!exists(dirname)) {
    if(!create_directory(dirname)) { std::cout << " error: couldn't create directory" << dirname << "\n";}
  }
}

void o(vector<vector<double>> vvd) {
  for(auto i: vvd) { boost::copy(i, ostream_iterator<double>(cout<<fixed," ")); cout <<"\n";}
}
void o(vector<vector<double>> vvd, string ofn) {
  ofstream of(ofn);
  for(auto i: vvd) { boost::copy(i, ostream_iterator<double>(of<<fixed," ")); of <<"\n";}
}

void o(vector<double> v, string ofn) {
  ofstream of(ofn);
  boost::copy(v,ostream_iterator<double>(of,"  ")); of << "\n";
}

void o(vector<double> v) {
  boost::copy(v,ostream_iterator<double>(cout,"  ")); cout << "\n";
}

void o(string name,vector<double> v) {
  cout << name << " -->  ";
  o(v);
}
