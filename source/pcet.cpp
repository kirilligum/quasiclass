#include <vector>
#include <algorithm>

#include "pcet.hpp"

using namespace std;


//pcet::pcet( vector<double> param, vector<double> wbath, vector<double> cbath): 
  //mp(param[0]),
  //w1p(param[1]),
  //w2p(param[2]),
  //r1p(param[3]),
  //r2p(param[4]),
  //e1(param[5]),
  //e2(param[6]),
  //delta(param[7]),
  //n_shift(param[8]),
  //1_div_2mp(0.5*mp),
  //wb(wbath), cb(cbath) {}

//double pcet::u(double qp, int i) {
  //if(i!=1 && i!=2) cout << "error \n";
  //return mp*wp[i]*wp[i]*0.5*(qp-rp[i])*(qp-rp[i]),
//}

//double pcet::vc() {
  //return 2.*delta*sqrt((n1+n_coup_shift)*(n2+n_coup_shift))*xxpp;
//}

//double pcet::energy(const vector<double> &s){
  //double t = s[0],
         //pp = s[1],
         //qp = s[2],
         //p1 = s[3],
         //x1 = s[4],
         //p2 = s[5],
         //x2 = s[6];
  //double tp = pp*pp*1_div_2mr,
         //u1 = u(qp,1);
         //u1 = u(qp,2);
         //n1 = get_n{n_shift}(s,1),
         //n2 = get_n{n_shift}(s,2),
         //v1 = (e1 + u1)*n1,
         //v2 = (e2 + u2)*n2,
         //xxpp = x1*x2 + p1*p2,
         //vc = 2.*delta*sqrt((n1+n_coup_shift)*(n2+n_coup_shift))*xxpp;
  //double hb = 0;
  //double vbcn=0.0;
  //double vbnn=0.0;
  //for (int i = 7; i < s.size()-1; ++(++i)) {
    //int ib=(i-7)/2;
    //vbcn +=2.0*cb[ib]*s[i+1]//check cb coef;
    //vbnn +=cb[ib]*cb[i];//check cb coef
    //tb+= s[i]*s[i];
    //vb+= s[i+1]*s[i+1]*wb[ib]*wb[ib];
  //}
  //tb*=0.5;
  //vb*=0.5;
  //vbcn*=n1;
  //vbnn*=n1*n1;
  //return tp+v1+v2+vc+tb+vb+vbcn+vbnn;
//}

//double pcet::du(double qp, int i) {
  //if(i!=1 && i!=2) cout << "error \n";
  //return mp*wp[i]*wp[i]*(qp-rp[i]),
//}

//double pcet::dn(char ch,int i) {
  //switch (ch) {
    //case 'x':
      //return x[i];
    //case 'p':
      //return p[i];
  //}
//}

//double pcet::dvn(char co, int i) {
  //return (e[i]+u[i])*dn(co,i);
//}

//double pcet::dxxpp(double ch, int i) {
  //switch (ch) {
    //case 'x':
      //return x[i+1];
    //case 'p':
      //return p[i+1];
  //}
//}

//double pcet::dvcn(char ch, int i) {
  //return -delta/sqrt((n(i)+n_coup_shift)*(n(i+1)+n_coup_shift))*(n(i+1)+n_coup_shift)*dn(ch,i)*xxpp +
  //2.*delta*sqrt((n(i)+n_coup_shift)*(n(i+1)+n_coup_shift))*dxxpp(ch,i);
//}

//double pcet::dvbcn(char ch, int i) {
  //return 0.0;
//}

//double pcet::dvbcn(char ch, int i) {
  //return 0.0;
//}

//void eom::operator() (const vector<double> &s, vector<double> &ds, const double){

  //ds[1] = -( du(qp,1)*n1+du(qp,2)*n2 ); // dp/dt = - dh/dq
  //ds[2] = pr*1_div_m; // dq/dt = dh/dp
  //ds[3] = -( dvn("x",1) +  dvcn("x",1) + dvbcn("x",1) + dvbnn("x",1)); // dp1/dt = -dh/dx1
  //ds[4] = dvn("p",1) +  dvcn("p",1) + dvbcn("p",1) + dvbnn("p",1); // dx1/dt = dh/dp1
  //ds[5] = -( dvn("x",2) +  dvcn("x",2)); // dp2/dt = -dh/dx2
  //ds[6] = dvn("p",2) +  dvcn("p",2); // dx2/dt = dh/dp2

  //for (int i = 5; i < s.size()-1; ++(++i)) {
    //int ib=(i-5)/2;
    //ds[i] = -( dvb_div_dq(i+1)+dvbcn_div_dq(i+1) );// dpb/dt = -dh/dqb
    //ds[i+1] = dtb_div_dp(i); // dqb/dt = dh/dpb
  //}
  ////sum_cq +=eps; // not sure about this yet. Steve has it

  //ds[0]=0.0;
  //ds.back()=0.0;

  //double nmn2 = 0.5*( s[1]*s[1] + s[2]*s[2] - s[3]*s[3] - s[4]*s[4]);
  //double sum_cq = 0.0;

  //for (int i = 5; i < s.size()-1; ++(++i)) {
    //int ib=(i-5)/2;
    //sum_cq += cb[ib]* s[i+1];
    //ds[i]=-(wb[ib]*wb[ib]*s[i+1]+cb[ib]*nmn2);
    //ds[i+1]= s[i];
  //}
  //sum_cq +=eps; // not sure about this yet. Steve has it

  //ds[1]=-(d*s[4]+s[2]*sum_cq);
  //ds[2]=(d*s[3]+s[1]*sum_cq);
  //ds[3]=-(d*s[2]-s[4]*sum_cq);
  //ds[4]=(d*s[1]-s[3]*sum_cq);


//}
