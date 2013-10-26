#pragma once
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
using namespace std;

inline double n_fromx(double p, double x) {
  return 0.5*(x*x+p*p-1.);
}

inline double q_fromx(double p, double x) {
  return atan2(p,x);
}


struct pcet {
  //using namespace std;
  const vector<double> bw,bc;
  const double m,e1,e2,w1,w2,rp1,rp2,v12,ns;
  pcet(const map<string,double> &system_parameters, const vector<double> &bath_w,  const vector<double> &bath_c) 
    : bw(bath_w),bc(bath_c),
      m       (system_parameters.at("m")),
      e1      (system_parameters.at("e1")),
      e2      (system_parameters.at("e2")),
      w1      (system_parameters.at("w1")),
      w2      (system_parameters.at("w2")),
      rp1     (system_parameters.at("rp1")),
      rp2     (system_parameters.at("rp2")),
      ns     (system_parameters.at("n_shift")),
      v12     (system_parameters.at("v12"))
  { }
  double operator() (const vector<double> &s) {
    double p = s[1];
    double r = s[2];
    double n1 = n_fromx(s[3],s[4]);
    double q1 = q_fromx(s[3],s[4]);
    double n2 = n_fromx(s[5],s[6]);
    double q2 = q_fromx(s[5],s[6]);
    double p1 = s[3];
    double x1 = s[4];
    double p2 = s[5];
    double x2 = s[6];
    //double pb = s[7];
    //double qb = s[8];

    //auto sys_ham = 
      //p*p*0.5/m + 
      //(e1+m*w1*w1*(r-rp1)*(r-rp1)*0.5)*0.5*(x1*x1+p1*p1-1.) + 
      //(e2+m*w2*w2*(r-rp2)*(r-rp2)*0.5)*0.5*(x2*x2+p2*p2-1.) +
      //2.0 * v12 * sqrt((0.5*(x1*x1+p1*p1-1.)+0.5)*(0.5*(x2*x2+p2*p2-1.)+0.5)) *
        //cos(atan2(p1,x1)-atan2(p2,x2));

    //double  sum_hamc =0;
    //for (size_t i = 0; i < bw.size(); ++i) {
      //auto pb = s[7+2*i];
      //auto qb = s[8+2*i];
      //sum_hamc+= 
        //pow(pb,2)/2. + 
        //(pow(bw[i],2)*
         //pow(qb + (bc[i]*(-ns + pow(p1,2) + pow(x1,2)))/(2.*pow(bw[i],2)),2)
         //)/2.;
    //}
    //auto ham = (0.5*pow(p,2))/m + 
      //0.5*(e1 + 0.5*m*pow(r - rp1,2)*pow(w1,2))*(-ns + pow(p1,2) + pow(x1,2)) + 
      //sum_hamc +
      //0.5*(e2 + 0.5*m*pow(r - rp2,2)*pow(w2,2))*(-ns + pow(p2,2) + pow(x2,2)) + 
      //2.*v12*(p1*p2 + x1*x2)*
        //sqrt(
            //(0.5 + 0.5*(-ns + pow(p1,2) + pow(x1,2)))*
            //(0.5 + 0.5*(-ns + pow(p2,2) + pow(x2,2))))
        //+
      //0.0
      ;
    //double p1 = s[1];
    //double x1 = s[2];
    //double p2 = s[3];
    //double x2 = s[4];
    
    auto nnd2 = 0.5*(x1*x1 + p1*p1 - x2*x2 - p2*p2);

    auto ecoup = v12*(x1*x2 + p1*p2);

    double bath_ham=0;
    double ebath_coup=0;
    for (size_t i = 0; i < bw.size(); ++i) {
      auto pb = s[7+2*i];
      auto qb = s[8+2*i];
      bath_ham+= pb*pb + bw[i]*bw[i]*qb*qb;
      ebath_coup+= bc[i]*qb;
    }
    bath_ham*=0.5;
    ebath_coup*=nnd2;

    //return (ecoup);
    auto ham = (bath_ham + ebath_coup + ecoup); // may be - ebath_coup

    //auto sys_ham = p*p*0.5 + r*r*0.5;
    //auto sys_ham = p*p*0.5/m + (e1+u1(r))*n1 + (e2+u2(r))*n2;
    //auto sys_ham = p*p*0.5/m + (e1+u1(r))*n1 + (e2+u2(r))*n2 +2.0*v12*sqrt((n1+0.5)*(n2+0.5))*cos(q1-q2);
    //auto bath_ham = pb*pb*0.5/mb + mb*wb*wb*0.5*(qb+n1*c/(mb*wb*wb))*(qb+n1*c/(mb*wb*wb));
    //double bath_ham;
    //for (size_t i = 0; i < bw.size(); ++i) {
      //auto pb = s[7+2*i];
      //auto qb = s[8+2*i];
      //bath_ham+= pb*pb*0.5 
        //+ bw[i]*bw[i]*0.5*(
            //qb+n1*bc[i]/(bw[i]*bw[i]))*(
            //qb+n1*bc[i]/(bw[i]*bw[i]));
    //}
    //return (sys_ham + bath_ham);
    return ham;
  }
  //double m, e1, e2, m, m, w1, w2, rp1, rp2, v12;
  double u1(double r) { 
    return m*w1*w1*(r-rp1)*(r-rp1)*0.5;}
  double u2(double r) { 
    return m*w2*w2*(r-rp2)*(r-rp2)*0.5;}
};

struct d_pcet : pcet {
  using pcet::pcet;
  vector<double> operator() (const vector<double> &s) ;
  double d_u1_dr(double r) { 
    return m*w1*w1*(r-rp1);
  }
  double d_u2_dr(double r) { 
    return m*w2*w2*(r-rp2);
  }
};

struct dfdx {
  function<double(double)> f_;
  double h_;
  dfdx(std::function<double(double)> f, double h) : f_(f), h_(h) {}
  double operator() (double x) {
    return (f_(x+h_)-f_(x-h_))/h_*0.5;
  }
};

inline vector<double> dfdv_n (function<double(vector<double>)> f, vector<double> &v, double h) { ///> can be modified to a structure and make an option with a copy of v
  vector<double> vp(v);
  boost::for_each(vp,[h](double x) { return x+=h;});
  vector<double> vm(v);
  boost::for_each(vm,[h](double x) { return x-=h;});
  std::vector<double> fv (v.size(),0.0);
  boost::for_each(v,fv,[&v,h,f](double &x, double &fi) {
      auto tmp = x;
      x+=h;
      auto fp =f(v);
      x = tmp -h;
      auto fm =f(v);
      x = tmp;
      fi =  (fp-fm)/h*0.5;
      });
  return fv;
}

inline double dfdv_n (function<double(vector<double>)> f, vector<double> &v, int i, double h) { ///> can be modified to a structure and make an option with a copy of v
  vector<double> vp(v);
  boost::for_each(vp,[h](double x) { return x+=h;});
  vector<double> vm(v);
  boost::for_each(vm,[h](double x) { return x-=h;});
  std::vector<double> fv (v.size(),0.0);
  boost::for_each(v,fv,[&v,h,f](double &x, double &fi) {
      auto tmp = x;
      x+=h;
      auto fp =f(v);
      x = tmp -h;
      auto fm =f(v);
      x = tmp;
      fi =  (fp-fm)/h*0.5;
      });
  return fv[i];
}

inline vector<double> dfdv_n (function<double(vector<double>)> f, const vector<double> &cv, double h) { ///> can be modified to a structure and make an option with a copy of v
  auto v = cv;
  return dfdv_n(f,v,h);
}

struct dfdv_num {
  function<double(vector<double>)> f;
  double h;
  dfdv_num (function<double(vector<double>)> f_, double h_) 
    : f(f_), h(h_) {}
  vector<double> operator()(const vector<double> &cv) {
    //auto v= cv;
    return dfdv_n(f,cv,h);
  }
};

struct eom {
  function<double(vector<double>)> f_;
  function<vector<double>(vector<double>)> dfdv;
  //double h_;
  eom( function<vector<double>(vector<double>)> dfdv_) : dfdv(dfdv_) {}
  //eom( function<double(vector<double>)> f, double h) : f_(f), h_(h) {}
  void operator () ( const vector<double> &v, vector<double> &dvdt, const double ) {
    auto dfdv_at_v = dfdv(v);
    //auto dfdv_at_v = dfdv(f_,v,h_);
    dvdt[0] = 0.0; ///> the first element is time that will be calculated later
    for (size_t i = 1; i < dfdv_at_v.size(); ++(++i)) {
      dvdt[i] = -dfdv_at_v[i+1];
      dvdt[i+1] = dfdv_at_v[i];
    }
  }
};

struct streamer {
  std::vector<std::vector<double>>& states_;
  streamer(std::vector<std::vector<double>>& states) : states_(states) {}
  void operator()(const std::vector<double>&x, double t) {
    states_.emplace_back(x);
    //states_.back().insert(states_.back().begin(),t);
    states_.back().front() = t;
  }
};


//struct pcet {
  //double mp, delta, n_shift, 1_div_2mp;
  //vector<double> param,wb,cb;
  //vector<double> wp,rp,e;
  //pcet(vector<double> param_, vector<double> wb_, vector<double> cb_);
  //double energy(const vector<double> &s);
  //void eom(const vector<double> &s, vector<double> &ds, const double);
//};
