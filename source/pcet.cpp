#include <vector>
#include <algorithm>

#include "pcet.hpp"
#include <vector>
#include <map>
//#include <tuple>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>


vector<double> d_pcet::operator() (const vector<double> &s) {
  using std::cout;
  using std::vector;
  using namespace boost;
  using namespace boost::adaptors;
  typedef vector<double> vd; typedef vector<vector<double>> vvd; typedef vector<vector<vector<double>>> vvvd;
  vector<double> ds(s.size());
  //double p = s[1];
  //double r = s[2];
  //double p1 = s[3];
  //double x1 = s[4];
  //double p2 = s[5];
  //double x2 = s[6];

  double nmn2 = 0.5*( s[3]*s[3] + s[4]*s[4] - s[5]*s[5] - s[6]*s[6]);
  double sum_cq = 0.0;
  double d=v12;
//for(auto & i :bc) {cout << i << "  " ;}
//cout << "\n";

  //cout << bc.size() << "b";
  //cout << s.size() << "s";
  for (int ib = 0; ib < bc.size(); ++ib) {
  //for (int i = 7; i < s.size(); ++(++i)) {
    //int ib=(i-7)/2;
    //double c=0.0;
    //cout << i<< "." << ib << "_";
    double c=bc[ib];
    int i = 2*ib+7;
    sum_cq += c* s[i+1];
    //sum_cq += c* s[i+1];
    //sum_cq += bc[ib]* s[i+1];
    //ds[i]=-(bw[ib]*bw[ib]*s[i+1]+bc[ib]*nmn2);
    //ds[i+1]= s[i];
    ds[i]= s[i];
    ds[i+1]=(bw[ib]*bw[ib]*s[i+1]+c*nmn2);
    //ds[i+1]=(bw[ib]*bw[ib]*s[i+1]+bc[ib]*nmn2);
  }
  //sum_cq +=eps; // not sure about this yet. Steve has it
  //cout << "{" << sum_cq << "}";

  ds[3]=(d*s[5]+s[3]*sum_cq);
  ds[4]=(d*s[6]+s[4]*sum_cq);
  ds[5]=(d*s[3]-s[5]*sum_cq);
  ds[6]=(d*s[4]-s[6]*sum_cq);

  //ds[3]=-(d*s[6]+s[4]*sum_cq);
  //ds[4]=(d*s[5]+s[3]*sum_cq);
  //ds[5]=-(d*s[4]-s[6]*sum_cq);
  //ds[6]=(d*s[3]-s[5]*sum_cq);

  ds[0]=0.0;
  //ds.back()=0.0;

  //double sum_dp1c = 0, sum_dx1c = 0, sum_hamc =0;
  //for (size_t i = 0; i < bw.size(); ++i) {
    //auto pb = s[7+2*i];
    //auto qb = s[8+2*i];
    //auto d_sys_ham_dpb = pb;
    //auto d_sys_ham_dqb = 1.*pow(bw[i],2)*(qb + 0.5*(bc[i]*(-ns + pow(p1,2) + pow(x1,2)))/pow(bw[i],2));
    //ds[7+2*i] = d_sys_ham_dpb;
    //ds[8+2*i] = d_sys_ham_dqb;
    //sum_dp1c+= bc[i]*p1*(qb + (bc[i]*(-ns + pow(p1,2) + pow(x1,2)))/(2.*pow(bw[i],2))) ;
    //sum_dx1c+= bc[i]*x1*(qb + (bc[i]*(-ns + pow(p1,2) + pow(x1,2)))/(2.*pow(bw[i],2))) ;
    //sum_hamc+= 
      //pow(pb,2)/2. + 
      //(pow(bw[i],2)*pow(qb + (bc[i]*(-ns + pow(p1,2) + pow(x1,2)))/(2.*pow(bw[i],2)),2))/2.;
  //}

  //auto ham = (0.5*pow(p,2))/m + 0.5*(e1 + 0.5*m*pow(r - rp1,2)*pow(w1,2))*(-ns + pow(p1,2) + pow(x1,2)) + 
    //sum_hamc +
 //0.5*(e2 + 0.5*m*pow(r - rp2,2)*pow(w2,2))*(-ns + pow(p2,2) + pow(x2,2)) + 2.*v12*(p1*p2 + x1*x2)*sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))));

  //auto nnd2 = 0.5*(x1*x1 + p1*p1 - x2*x2 - p2*p2);

  //auto ecoup = v12*(x1*x2 + p1*p2);

  //double bath_ham=0;
  //double ebath_coup=0;
  //for (size_t i = 0; i < bw.size(); ++i) {
    //auto pb = s[7+2*i];
    //auto qb = s[8+2*i];
    //bath_ham+= pb*pb + bw[i]*bw[i]*qb*qb;
    //ebath_coup+= bc[i]*qb;
    //auto d_sys_ham_dpb = pb;
    //auto d_sys_ham_dqb = pow(bw[i],2)*qb;
    //ds[7+2*i] = d_sys_ham_dpb;
    //ds[8+2*i] = d_sys_ham_dqb;
  //}
  //bath_ham*=0.5;
  //ebath_coup*=nnd2;


  //return (ecoup);
  //auto ham = (bath_ham - ebath_coup + ecoup);
  //auto ham = (bath_ham + ebath_coup + ecoup); // may be - ebath_coup

  //auto d_sys_ham_dp = 0.0;
  //auto d_sys_ham_dr = 0.0;
  //auto d_sys_ham_dp1 =  p1*ebath_coup + v12*p2 ;
  //auto d_sys_ham_dx1 =  x1*ebath_coup + v12*x2 ;
  //auto d_sys_ham_dp2 = -p2*ebath_coup + v12*p1 ;
  //auto d_sys_ham_dx2 = -x2*ebath_coup + v12*x1 ;

  //auto d_sys_ham_dp = (1.*p)/m;
  //auto d_sys_ham_dr = 
    //0.5*m*(r - rp1)*pow(w1,2)*(-ns + pow(p1,2) + pow(x1,2)) + 
    //0.5*m*(r - rp2)*pow(w2,2)*(-ns + pow(p2,2) + pow(x2,2));
  //auto d_sys_ham_dp1 = 0.0;
  //auto d_sys_ham_dx1 = 0.0;
  //auto d_sys_ham_dp2 = 0.0;
  //auto d_sys_ham_dx2 = 0.0;

  //auto d_sys_ham_dp1 = 1.*p1*(e1 + 0.5*m*pow(r - rp1,2)*pow(w1,2)) + 
    //sum_dp1c +
 //(1.*p1*v12*(p1*p2 + x1*x2)*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))))/sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2)))) + 
 //2.*p2*v12*sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))));
  //auto d_sys_ham_dx1 = 1.*(e1 + 0.5*m*pow(r - rp1,2)*pow(w1,2))*x1 + 
    //sum_dx1c +
 //(1.*v12*x1*(p1*p2 + x1*x2)*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))))/sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2)))) + 
 //2.*v12*x2*sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))));
  //auto d_sys_ham_dp2 = 1.*p2*(e2 + 0.5*m*pow(r - rp2,2)*pow(w2,2)) + (1.*p2*v12*(0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(p1*p2 + x1*x2))/
  //sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2)))) + 2.*p1*v12*sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))));
  //auto d_sys_ham_dx2 = 1.*(e2 + 0.5*m*pow(r - rp2,2)*pow(w2,2))*x2 + (1.*v12*(0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*x2*(p1*p2 + x1*x2))/
  //sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2)))) + 2.*v12*x1*sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))));

  //auto d_sys_ham_dp = (1.*p)/m;
  //auto d_sys_ham_dr = 0.5*m*(r - rp1)*pow(w1,2)*(-ns + pow(p1,2) + pow(x1,2)) + 0.5*m*(r - rp2)*pow(w2,2)*(-ns + pow(p2,2) + pow(x2,2));
  //auto d_sys_ham_dp1 = 1.*p1*(e1 + 0.5*m*pow(r - rp1,2)*pow(w1,2)) + 
    //sum_dp1c +
 //(1.*p1*v12*(p1*p2 + x1*x2)*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))))/sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2)))) + 
 //2.*p2*v12*sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))));
  //auto d_sys_ham_dx1 = 1.*(e1 + 0.5*m*pow(r - rp1,2)*pow(w1,2))*x1 + 
    //sum_dx1c +
 //(1.*v12*x1*(p1*p2 + x1*x2)*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))))/sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2)))) + 
 //2.*v12*x2*sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))));
  //auto d_sys_ham_dp2 = 1.*p2*(e2 + 0.5*m*pow(r - rp2,2)*pow(w2,2)) + (1.*p2*v12*(0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(p1*p2 + x1*x2))/
  //sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2)))) + 2.*p1*v12*sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))));
  //auto d_sys_ham_dx2 = 1.*(e2 + 0.5*m*pow(r - rp2,2)*pow(w2,2))*x2 + (1.*v12*(0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*x2*(p1*p2 + x1*x2))/
  //sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2)))) + 2.*v12*x1*sqrt((0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))));

  //auto sys_ham = 
    //p*p*0.5/m 
    //+ 
    //(e1+m*w1*w1*(r-rp1)*(r-rp1)*0.5)*0.5*(x1*x1+p1*p1-1.) 
    //+ 
    //(e2+m*w2*w2*(r-rp2)*(r-rp2)*0.5)*0.5*(x2*x2+p2*p2-1.) 
    //+
    //2.0 * v12 * sqrt((0.5*(x1*x1+p1*p1-1.)+0.5)*(0.5*(x2*x2+p2*p2-1.)+0.5)) 
      //*
      //cos(atan2(p1,x1)-atan2(p2,x2))
    //;

  //auto d_sys_ham_dp = p/m ;
  //ds[1]= d_sys_ham_dp;

  //auto d_sys_ham_dr =
    //(m*w1*w1*(r-rp1))*0.5*(x1*x1+p1*p1-1.) + 
    //(m*w2*w2*(r-rp2))*0.5*(x2*x2+p2*p2-1.)
    //;
  //ds[2]= d_sys_ham_dr;

  //auto d_sys_ham_dp1 = 
    //p1*(e1 + 0.5*m*pow(r - rp1,2)*pow(w1,2)) + 
     //(1.*p1*v12*(0.5 + 0.5*(-1. + pow(p2,2) + pow(x2,2)))*
        //cos(atan2(p1,x1) - atan2(p2,x2)))/
      //sqrt((0.5 + 0.5*(-1. + pow(p1,2) + pow(x1,2)))*
        //(0.5 + 0.5*(-1. + pow(p2,2) + pow(x2,2)))) + 
     //(2.*v12*x1*sqrt((0.5 + 0.5*(-1. + pow(p1,2) + pow(x1,2)))*
          //(0.5 + 0.5*(-1. + pow(p2,2) + pow(x2,2))))*
        //sin(atan2(p1,x1) - atan2(p2,x2)))/
      //(pow(p1,2) + pow(x1,2));
  //for (size_t i = 0; i < bw.size(); ++i) {
    //auto pb = s[7+2*i];
    //auto qb = s[8+2*i];
    //d_sys_ham_dp1 += 1.*bc[i]*p1*(qb + (0.5*bc[i]*(-1 + pow(p1,2) + pow(x1,2)))/pow(bw[i],2));
  //}
  //ds[3] = d_sys_ham_dp1;

  //auto d_sys_ham_dx1 = 
    //1.*(e1 + 0.5*m*pow(r - rp1,2)*pow(w1,2))*x1 + 
     //(1.*v12*x1*(0.5 + 0.5*(-1. + pow(p2,2) + pow(x2,2)))*
        //cos(atan2(p1,x1) - atan2(p2,x2)))/
      //sqrt((0.5 + 0.5*(-1. + pow(p1,2) + pow(x1,2)))*
        //(0.5 + 0.5*(-1. + pow(p2,2) + pow(x2,2)))) - 
     //(2.*p1*v12*sqrt((0.5 + 0.5*(-1. + pow(p1,2) + pow(x1,2)))*
          //(0.5 + 0.5*(-1. + pow(p2,2) + pow(x2,2))))*
        //sin(atan2(p1,x1) - atan2(p2,x2)))/
      //(pow(p1,2) + pow(x1,2));
  //for (size_t i = 0; i < bw.size(); ++i) {
    //auto pb = s[7+2*i];
    //auto qb = s[8+2*i];
    //d_sys_ham_dx1 += 1.*bc[i]*x1*(qb + (0.5*bc[i]*(-1 + pow(p1,2) + pow(x1,2)))/pow(bw[i],2));
  //}
  //ds[4] = d_sys_ham_dx1;

  //auto d_sys_ham_dp2 = 
    //p2*(e2 + 0.5*m*pow(r - rp2,2)*pow(w2,2)) + 
     //(1.*p2*v12*(0.5 + 0.5*(-1. + pow(p1,2) + pow(x1,2)))*
        //cos(atan2(p1,x1) - atan2(p2,x2)))/
      //sqrt((0.5 + 0.5*(-1. + pow(p1,2) + pow(x1,2)))*
        //(0.5 + 0.5*(-1. + pow(p2,2) + pow(x2,2)))) - 
     //(2.*v12*x2*sqrt((0.5 + 0.5*(-1. + pow(p1,2) + pow(x1,2)))*
          //(0.5 + 0.5*(-1. + pow(p2,2) + pow(x2,2))))*
        //sin(atan2(p1,x1) - atan2(p2,x2)))/
      //(pow(p2,2) + pow(x2,2));
  //ds[5] = d_sys_ham_dp2;

  //auto d_sys_ham_dx2 = 
    //1.*(e2 + 0.5*m*pow(r - rp2,2)*pow(w2,2))*x2 + 
     //(1.*v12*(0.5 + 0.5*(-1. + pow(p1,2) + pow(x1,2)))*x2*cos(atan2(p1,x1) - atan2(p2,x2)))/
      //sqrt((0.5 + 0.5*(-1. + pow(p1,2) + pow(x1,2)))*
        //(0.5 + 0.5*(-1. + pow(p2,2) + pow(x2,2)))) + 
     //(2.*p2*v12*sqrt((0.5 + 0.5*(-1. + pow(p1,2) + pow(x1,2)))*
          //(0.5 + 0.5*(-1. + pow(p2,2) + pow(x2,2))))*sin(atan2(p1,x1) - atan2(p2,x2)))/
      //(pow(p2,2) + pow(x2,2));
  //ds[6] = d_sys_ham_dx2;

  //for (size_t i = 0; i < bw.size(); ++i) {
    //auto pb = s[7+2*i];
    //auto qb = s[8+2*i];
    //auto d_sys_ham_dpb = pb;
    //auto d_sys_ham_dqb = 1.*pow(bw[i],2)*(qb + (0.5*bc[i]*(-1 + pow(p1,2) + pow(x1,2)))/pow(bw[i],2));
    //ds[7+2*i] = d_sys_ham_dpb;
    //ds[8+2*i] = d_sys_ham_dqb;
  //}

  return ds;
}
