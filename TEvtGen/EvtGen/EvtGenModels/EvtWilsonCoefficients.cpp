//---------------------------------------------------------------------------------
//
// Wilson coeficients according to A.J.Buras and M.Munz, Phys.Rev. D52, 186. (1995)
// Thanks to N. Nikitine for example code for Pythia
// Coefficient C8eff and C2 correction to C7eff taken from:
//   A.J.Buras, M.Misiak, M.Munz, S.Pokorski, Nucl.Phys. B424, 374 (1994)
//
// Used constants come from PDG 2004
//
// P. Reznicek 18.02.2005
//
//  04/03/2005  PR   Added h-function
//
//---------------------------------------------------------------------------------

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtWilsonCoefficients.hh"
#include <stdlib.h>

// including EvtLi2Spence.F
extern "C" {
  extern double li2spence_(double*);
}

EvtWilsonCoefficients::EvtWilsonCoefficients(){
  int i,j;
  double tmpa[8]={14./23.,16./23.,6./23.,-12./23.,0.4086,-0.4230,-0.8994,0.1456};
  double tmph[8]={2.2996,-1.0880,-3./7.,-1./14.,-0.6494,-0.0380,-0.0186,-0.0057};
  double tmpp[8]={0,0,-80./203.,8./33.,0.0433,0.1384,0.1648,-0.0073};
  double tmps[8]={0,0,-0.2009,-0.3579,0.0490,-0.3616,-0.3554,0.0072};
  double tmpq[8]={0,0,0,0,0.0318,0.0918,-0.2700,0.0059};
  double tmpg[8]={313063./363036.,0,0,0,-0.9135,0.0873,-0.0571,-0.0209};
  double tmpk[6][8]={{0,0,1./2.,-1./2.,0,0,0,0},
                     {0,0,1./2.,+1./2.,0,0,0,0},
	             {0,0,-1./14.,+1./6.,0.0510,-0.1403,-0.0113,0.0054},
  	             {0,0,-1./14.,-1./6.,0.0984,+0.1214,+0.0156,0.0026},
	             {0,0,0,0,-0.0397,0.0117,-0.0025,+0.0304},
	             {0,0,0,0,+0.0335,0.0239,-0.0462,-0.0112}};
  double tmpr[2][8]={{0,0,+0.8966,-0.1960,-0.2011,0.1328,-0.0292,-0.1858},
                     {0,0,-0.1193,+0.1003,-0.0473,0.2323,-0.0133,-0.1799}};
  for(i=0;i<8;i++){
    a[i]=tmpa[i];
    h[i]=tmph[i];
    p[i]=tmpp[i];
    s[i]=tmps[i];
    q[i]=tmpq[i];
    g[i]=tmpg[i];
    for(j=0;j<6;j++) k[j][i]=tmpk[j][i];
    for(j=0;j<2;j++) r[j][i]=tmpr[j][i];
  }
  m_n_f=5;
  m_Lambda=0.2167;
  m_alphaMZ=0.1187;
  m_mu=4.8;
  m_M_Z=91.1876;
  m_M_t=174.3;
  m_M_W=80.425;
  m_alphaS=0;
  m_eta=0;
  m_sin2W=0.23120;
  m_ialpha=137.036;
  m_C1=m_C2=m_C3=m_C4=m_C5=m_C6=m_C7=m_C7eff0=m_C8=m_C8eff0=m_C9=m_C9tilda=m_C10=m_C10tilda=m_P0=0;
  m_A=m_B=m_C=m_D=m_E=m_F=m_Y=m_Z=m_PE=0;
  m_ksi=0;
}

void EvtWilsonCoefficients::SetRenormalizationScheme(std::string scheme){
  if(scheme=="NDR") m_ksi=0;
  else if(scheme=="HV") m_ksi=1;
  else{
    report(Severity::Error,"EvtGen") << "ERROR: EvtWilsonCoefficients knows only NDR and HV schemes !" << std::endl;
    ::abort();
  }
}

double EvtWilsonCoefficients::alphaS(double mu=4.8,int n_f=5,double Lambda=0.2167){
// calculate strong coupling constant for n_f flavours and scale mu
  double beta0=11.-2./3.*n_f;
  double beta1=51.-19./3.*n_f;
  double beta2=2857.-5033./9.*n_f+325./27.*n_f*n_f;
  double lnratio=log(mu*mu/Lambda/Lambda);
  double aS=4.*EvtConst::pi/beta0/lnratio*(1.-2*beta1/beta0/beta0*log(lnratio)/lnratio+
            4*beta1*beta1/beta0/beta0/beta0/beta0/lnratio/lnratio*((log(lnratio)-0.5)*(log(lnratio)-0.5)+beta2*beta0/8/beta1/beta1-5./4.));
  return aS;
}

double EvtWilsonCoefficients::Lambda(double alpha=0.1187,int n_f=5,double mu=91.1876,double epsilon=0.00005,int maxstep=1000){
// calculate Lambda matching alphaS using simple iterative method
  int i;
  double difference=0;
  double Lambda=mu*0.9999999999;
  double step=-mu/20;
  for(i=0;i<maxstep && (difference=fabs(alphaS(mu,n_f,Lambda)-alpha))>=epsilon;i++){
    report(Severity::Debug,"EvtGen") << " Difference of alpha_S from " << alpha << " is " << difference << " at Lambda = " << Lambda << std::endl;
    if(alphaS(mu,n_f,Lambda)>alpha){
      if(step>0) step*=-0.4;
      if(alphaS(mu,n_f,Lambda+step-epsilon)<alphaS(mu,n_f,Lambda+step)) Lambda+=step;
      else step*=0.4;
    }else{
      if(step<0) step*=-0.4;
      if(Lambda+step<mu) Lambda+=step;
      else step*=0.4;
    }
  }
  report(Severity::Debug,"EvtGen") << " Difference of alpha_S from " << alpha << " is " << difference << " at Lambda = " << Lambda << std::endl;
  if(difference>=epsilon){
    report(Severity::Error,"EvtGen") << " ERROR: Did not converge Lambda for alpha_s = " << alpha << " , difference " << difference << " >= " << epsilon << " after " << i << " steps !" << std::endl;
    ::abort();
    return -1;
  }else{
    report(Severity::Info,"EvtGen") << " For alpha_s = " << alphaS(mu,n_f,Lambda) << " was found Lambda = " << Lambda << std::endl;
    return Lambda;
  }
}

double EvtWilsonCoefficients::eta(double mu=4.8,int n_f=5,double Lambda=0.2167,double M_W=80.425){
  return alphaS(M_W,n_f,Lambda)/alphaS(mu,n_f,Lambda);
}


EvtComplex EvtWilsonCoefficients::C1(double mu=4.8,int n_f=5,double Lambda=0.2167,double M_W=80.425){
  int i;
  EvtComplex myC1(0,0);
  for(i=0;i<8;i++) myC1+=k[0][i]*pow(eta(mu,n_f,Lambda,M_W),a[i]);
  return myC1;
}

EvtComplex EvtWilsonCoefficients::C2(double mu=4.8,int n_f=5,double Lambda=0.2167,double M_W=80.425){
  int i;
  EvtComplex myC2(0,0);
  for(i=0;i<8;i++) myC2+=k[1][i]*pow(eta(mu,n_f,Lambda,M_W),a[i]);
  return myC2;
}

EvtComplex EvtWilsonCoefficients::C3(double mu=4.8,int n_f=5,double Lambda=0.2167,double M_W=80.425){
  int i;
  EvtComplex myC3(0,0);
  for(i=0;i<8;i++) myC3+=k[2][i]*pow(eta(mu,n_f,Lambda,M_W),a[i]);
  return myC3;
}

EvtComplex EvtWilsonCoefficients::C4(double mu=4.8,int n_f=5,double Lambda=0.2167,double M_W=80.425){
  int i;
  EvtComplex myC4(0,0);
  for(i=0;i<8;i++) myC4+=k[3][i]*pow(eta(mu,n_f,Lambda,M_W),a[i]);
  return myC4;
}

EvtComplex EvtWilsonCoefficients::C5(double mu=4.8,int n_f=5,double Lambda=0.2167,double M_W=80.425){
  int i;
  EvtComplex myC5(0,0);
  for(i=0;i<8;i++) myC5+=k[4][i]*pow(eta(mu,n_f,Lambda,M_W),a[i]);
  return myC5;
}

EvtComplex EvtWilsonCoefficients::C6(double mu=4.8,int n_f=5,double Lambda=0.2167,double M_W=80.425){
  int i;
  EvtComplex myC6(0,0);
  for(i=0;i<8;i++) myC6+=k[5][i]*pow(eta(mu,n_f,Lambda,M_W),a[i]);
  return myC6;
}


EvtComplex EvtWilsonCoefficients::C7(double M_t=174.3,double M_W=80.425){
  return EvtComplex(-0.5*A(M_t*M_t/M_W/M_W),0);
}

EvtComplex EvtWilsonCoefficients::C8(double M_t=174.3,double M_W=80.425){
  return EvtComplex(-0.5*F(M_t*M_t/M_W/M_W),0);
}

EvtComplex EvtWilsonCoefficients::C7eff0(double mu=4.8,int n_f=5,double Lambda=0.2167,double M_t=174.3,double M_W=80.425){
  int i;
  EvtComplex myC7eff(0,0);
  for(i=0;i<8;i++) myC7eff+=h[i]*pow(eta(mu,n_f,Lambda,M_W),a[i]);
  myC7eff*=C2(mu,n_f,Lambda,M_W);
  myC7eff+=pow(eta(mu,n_f,Lambda,M_W),16./23.)*C7(M_t,M_W);
  myC7eff+=8./3.*(pow(eta(mu,n_f,Lambda,M_W),14./23.)-pow(eta(mu,n_f,Lambda,M_W),16./23.))*C8(M_t,M_W);
  return myC7eff;
}

EvtComplex EvtWilsonCoefficients::C8eff0(double mu=4.8,int n_f=5,double Lambda=0.2167,double M_t=174.3,double M_W=80.425){
  int i;
  EvtComplex myC8eff(0,0);
  for(i=0;i<8;i++) myC8eff+=g[i]*pow(eta(mu,n_f,Lambda,M_W),a[i]);
  myC8eff+=pow(eta(mu,n_f,Lambda,M_W),14./23.)*C8(M_t,M_W);
  return myC8eff;
}


EvtComplex EvtWilsonCoefficients::C10tilda(double sin2W=0.23120,double M_t=174.3,double M_W=80.425){
  return EvtComplex(-Y(M_t*M_t/M_W/M_W)/sin2W,0);
}

EvtComplex EvtWilsonCoefficients::C10(double sin2W=0.23120,double M_t=174.3,double M_W=80.425,double ialpha=137.036){
  return ( 1./2/EvtConst::pi/ialpha*C10tilda(sin2W,M_t,M_W) );
}


double EvtWilsonCoefficients::A(double x){
  return ( x*(8*x*x+5*x-7)/12/pow(x-1,3) + x*x*(2-3*x)*log(x)/2/pow(x-1,4) );
}

double EvtWilsonCoefficients::B(double x){
  return ( x/4/(1-x) + x/4/(x-1)/(x-1)*log(x) );
}

double EvtWilsonCoefficients::C(double x){
  return ( x*(x-6)/8/(x-1) + x*(3*x+2)/8/(x-1)/(x-1)*log(x) );
}

double EvtWilsonCoefficients::D(double x){
  return ( (-19*x*x*x+25*x*x)/36/pow(x-1,3) + x*x*(5*x*x-2*x-6)/18/pow(x-1,4)*log(x) - 4./9*log(x) );
}

double EvtWilsonCoefficients::E(double x){
  return ( x*(18-11*x-x*x)/12/pow(1-x,3) + x*x*(15-16*x+4*x*x)/6/pow(1-x,4)*log(x) - 2./3*log(x) );
}

double EvtWilsonCoefficients::F(double x){
  return ( x*(x*x-5*x-2)/4/pow(x-1,3) + 3*x*x/2/pow(x-1,4)*log(x) );
}

double EvtWilsonCoefficients::Y(double x){
  return (C(x)-B(x));
}

double EvtWilsonCoefficients::Z(double x){
  return (C(x)+1./4*D(x));
}


EvtComplex EvtWilsonCoefficients::C9(int ksi=0,double mu=4.8,int n_f=5,double Lambda=0.2167,double sin2W=0.23120,double M_t=174.3,double M_W=80.425,double ialpha=137.036){
  return ( 1./2/EvtConst::pi/ialpha*C9tilda(ksi,mu,n_f,Lambda,sin2W,M_t,M_W) );
}

EvtComplex EvtWilsonCoefficients::C9tilda(int ksi=0,double mu=4.8,int n_f=5,double Lambda=0.2167,double sin2W=0.23120,double M_t=174.3,double M_W=80.425){
  return ( P0(ksi,mu,n_f,Lambda,M_W) + Y(M_t*M_t/M_W/M_W)/sin2W - 4*Z(M_t*M_t/M_W/M_W) + PE(mu,n_f,Lambda,M_W)*E(M_t*M_t/M_W/M_W) );
}

EvtComplex EvtWilsonCoefficients::P0(int ksi=0,double mu=4.8,int n_f=5,double Lambda=0.2167,double M_W=80.425){
  int i;
  EvtComplex myP0(0,0);
  for(i=0;i<8;i++) myP0+=p[i]*pow(eta(mu,n_f,Lambda,M_W),a[i]+1);
  myP0=EvtConst::pi/alphaS(M_W,n_f,Lambda)*(-0.1875+myP0);
  myP0+=1.2468-ksi*4./9.*(3*C1(mu,n_f,Lambda,M_W)+C2(mu,n_f,Lambda,M_W)-C3(mu,n_f,Lambda,M_W)-3*C4(mu,n_f,Lambda,M_W));
  for(i=0;i<8;i++) myP0+=pow(eta(mu,n_f,Lambda,M_W),a[i])*(r[ksi][i]+s[i]*eta(mu,n_f,Lambda,M_W));
  return myP0;
}

double EvtWilsonCoefficients::PE(double mu=4.8,int n_f=5,double Lambda=0.2167,double M_W=80.425){
  int i;
  double myPE=0.1405;
  for(i=0;i<8;i++) myPE+=q[i]*pow(eta(mu,n_f,Lambda,M_W),a[i]+1);
  return myPE;
}


void EvtWilsonCoefficients::CalculateAllCoefficients(){
  m_Lambda=Lambda(m_alphaMZ,m_n_f,m_M_Z);
  m_C1=C1(m_mu,m_n_f,m_Lambda,m_M_W);
  m_C2=C2(m_mu,m_n_f,m_Lambda,m_M_W);
  m_C3=C3(m_mu,m_n_f,m_Lambda,m_M_W);
  m_C4=C4(m_mu,m_n_f,m_Lambda,m_M_W);
  m_C5=C5(m_mu,m_n_f,m_Lambda,m_M_W);
  m_C6=C6(m_mu,m_n_f,m_Lambda,m_M_W);
  m_C7=C7(m_M_t,m_M_W);
  m_C8=C8(m_M_t,m_M_W);
  m_C7eff0=C7eff0(m_mu,m_n_f,m_Lambda,m_M_t,m_M_W);
  m_C8eff0=C8eff0(m_mu,m_n_f,m_Lambda,m_M_t,m_M_W);
  m_C10tilda=C10tilda(m_sin2W,m_M_t,m_M_W);
  m_C10=C10(m_sin2W,m_M_t,m_M_W,m_ialpha);
  m_A=A(m_M_t*m_M_t/m_M_W/m_M_W);
  m_B=B(m_M_t*m_M_t/m_M_W/m_M_W);
  m_C=C(m_M_t*m_M_t/m_M_W/m_M_W);
  m_D=D(m_M_t*m_M_t/m_M_W/m_M_W);
  m_E=E(m_M_t*m_M_t/m_M_W/m_M_W);
  m_F=F(m_M_t*m_M_t/m_M_W/m_M_W);
  m_Y=Y(m_M_t*m_M_t/m_M_W/m_M_W);
  m_Z=Z(m_M_t*m_M_t/m_M_W/m_M_W);
  m_C9=C9(m_ksi,m_mu,m_n_f,m_Lambda,m_sin2W,m_M_t,m_M_W,m_ialpha);
  m_C9tilda=C9tilda(m_ksi,m_mu,m_n_f,m_Lambda,m_sin2W,m_M_t,m_M_W);
  m_P0=P0(m_ksi,m_mu,m_n_f,m_Lambda,m_M_W);
  m_PE=PE(m_mu,m_n_f,m_Lambda,m_M_W);
  m_alphaS=alphaS(m_mu,m_n_f,m_Lambda);
  m_eta=eta(m_mu,m_n_f,m_Lambda,m_M_W);
  report(Severity::Info,"EvtGen") << " +---------------------------------------" << std::endl;
  report(Severity::Info,"EvtGen") << " | Table of Wilson coeficients:" << std::endl;
  report(Severity::Info,"EvtGen") << " +---------------------------------------" << std::endl;
  report(Severity::Info,"EvtGen") << " | C1     =  " << m_C1 << std::endl;
  report(Severity::Info,"EvtGen") << " | C2     =  " << m_C2 << std::endl;
  report(Severity::Info,"EvtGen") << " | C3     =  " << m_C3 << std::endl;
  report(Severity::Info,"EvtGen") << " | C4     =  " << m_C4 << std::endl;
  report(Severity::Info,"EvtGen") << " | C5     =  " << m_C5 << std::endl;
  report(Severity::Info,"EvtGen") << " | C6     =  " << m_C6 << std::endl;
  report(Severity::Info,"EvtGen") << " | C7     =  " << m_C7 << std::endl;
  report(Severity::Info,"EvtGen") << " | C7eff0 =  " << m_C7eff0 << std::endl;
  report(Severity::Info,"EvtGen") << " | C8     =  " << m_C8 << std::endl;
  report(Severity::Info,"EvtGen") << " | C8eff0 =  " << m_C8eff0 << std::endl;
  report(Severity::Info,"EvtGen") << " | C9     =  " << m_C9 << std::endl;
  report(Severity::Info,"EvtGen") << " | C10    =  " << m_C10 << std::endl;
  report(Severity::Info,"EvtGen") << " +---------------------------------------" << std::endl;
  report(Severity::Info,"EvtGen") << " | Other constants:" << std::endl;
  report(Severity::Info,"EvtGen") << " +---------------------------------------" << std::endl;
  report(Severity::Info,"EvtGen") << " | Scale = " << m_mu << " GeV" << std::endl;
  report(Severity::Info,"EvtGen") << " | Number of effective flavors = " << m_n_f << std::endl;
  report(Severity::Info,"EvtGen") << " | Corresponding to aS(M_Z)" << "=" << m_alphaMZ << " Lambda = " << m_Lambda << " GeV" << std::endl;
  report(Severity::Info,"EvtGen") << " | Strong coupling constant = " << m_alphaS << std::endl;
  report(Severity::Info,"EvtGen") << " | Electromagnetic constant = 1/" << m_ialpha << std::endl;
  report(Severity::Info,"EvtGen") << " | Top mass = " << m_M_t << " GeV" << std::endl;
  report(Severity::Info,"EvtGen") << " | W-boson mass = " << m_M_W << " GeV" << std::endl;
  report(Severity::Info,"EvtGen") << " | Z-boson mass = " << m_M_Z << " GeV" << std::endl;
  report(Severity::Info,"EvtGen") << " | Sinus squared of Weinberg angle = " << m_sin2W << std::endl;
  report(Severity::Info,"EvtGen") << " +---------------------------------------" << std::endl;
  report(Severity::Debug,"EvtGen") << " | Intermediate functions:" << std::endl;
  report(Severity::Debug,"EvtGen") << " +---------------------------------------" << std::endl;
  report(Severity::Debug,"EvtGen") << " | A    = " << m_A << std::endl;
  report(Severity::Debug,"EvtGen") << " | B    = " << m_B << std::endl;
  report(Severity::Debug,"EvtGen") << " | C    = " << m_C << std::endl;
  report(Severity::Debug,"EvtGen") << " | D    = " << m_D << std::endl;
  report(Severity::Debug,"EvtGen") << " | E    = " << m_E << std::endl;
  report(Severity::Debug,"EvtGen") << " | F    = " << m_F << std::endl;
  report(Severity::Debug,"EvtGen") << " | Y    = " << m_Y << std::endl;
  report(Severity::Debug,"EvtGen") << " | Z    = " << m_Z << std::endl;
  report(Severity::Debug,"EvtGen") << " | eta  = " << m_eta << std::endl;
  report(Severity::Debug,"EvtGen") << " | C9~  = " << m_C9tilda << std::endl;
  report(Severity::Debug,"EvtGen") << " | C10~ = " << m_C10tilda << std::endl;
  report(Severity::Debug,"EvtGen") << " | P0   = " << m_P0 << std::endl;
  report(Severity::Debug,"EvtGen") << " | PE   = " << m_PE << std::endl;
  report(Severity::Debug,"EvtGen") << " +--------------------------------------" << std::endl;
}

EvtComplex EvtWilsonCoefficients::hzs(double z,double shat,double mu=4.8,double M_b=4.8){
  EvtComplex i1(0,1);
  double x=4.*z*z/shat;
  if(x==0)     return (8./27. - 8./9.*log(M_b/mu) - 4./9.*log(shat) + 4./9.*i1*EvtConst::pi);
  else if(x>1) return (8./27. - 8./9.*log(M_b/mu) - 8./9.*log(z) + 4./9.*x - 2./9.*(2.+x)*sqrt(x-1.) * 2*atan(1./sqrt(x-1.)));
  else         return (8./27. - 8./9.*log(M_b/mu) - 8./9.*log(z) + 4./9.*x - 2./9.*(2.+x)*sqrt(1.-x) * (log(fabs(sqrt(1.-x)+1)/fabs(sqrt(1.-x)-1))-i1*EvtConst::pi));
}

double EvtWilsonCoefficients::fz(double z){
  return (1. - 8.*z*z + 8.*pow(z,6.) - pow(z,8.) - 24.*pow(z,4.)*log(z));
}

double EvtWilsonCoefficients::kappa(double z,double alpha_S){
  return (1. - 2.*alpha_S/3./EvtConst::pi*((EvtConst::pi*EvtConst::pi-31./4.)*(1.-z)*(1.-z) + 1.5) );
}

double EvtWilsonCoefficients::etatilda(double shat,double alpha_S){
  return (1. + alpha_S/EvtConst::pi*omega(shat));
}

double EvtWilsonCoefficients::omega(double shat){
  double o=0;
  o -= (2./9.)*EvtConst::pi*EvtConst::pi;
  o -= (4./3.)*li2spence_(&shat);
  o -= (2./3.)*log(shat)*log(1.-shat);
  o -= log(1.-shat)*(5.+4.*shat)/(3.+6.*shat);
  o -= log(shat)*2.*shat*(1.+shat)*(1.-2.*shat)/3./(1.-shat)/(1.-shat)/(1.+2.*shat);
  o += (5.+9.*shat-6.*shat*shat)/6./(1.-shat)/(1.+2.*shat);
  return o;
}

EvtComplex EvtWilsonCoefficients::C9efftilda(double z,double shat,double alpha_S,EvtComplex c1,EvtComplex c2,EvtComplex c3,EvtComplex c4,EvtComplex c5,EvtComplex c6,EvtComplex c9tilda,int ksi=0){
  EvtComplex c(0,0);
  c += (c9tilda+ksi*4./9.*(3.*c1+c2-c3-3.*c4))*etatilda(shat,alpha_S);
  c += hzs(z,shat)*(3.*c1+c2+3.*c3+c4+3.*c5+c6);
  c -= 0.5*hzs(1,shat)*(4.*c3+4.*c4+3.*c5+c6);
  c -= 0.5*hzs(0,shat)*(c3+3.*c4);
  c += 2./9.*(3.*c3+c4+3.*c5+c6);
  return c;
}

EvtComplex EvtWilsonCoefficients::C7b2sg(double alpha_S,double et,EvtComplex c2,double M_t=174.3,double M_W=80.425){
  EvtComplex i1(0,1);
  return (i1*alpha_S*(2./9.*pow(et,14./23.)*(0.5*F(M_t*M_t/M_W/M_W)-0.1687)-0.03*c2));
}

EvtComplex EvtWilsonCoefficients::Yld(double q2,double *ki,double *Gi,double *Mi,int ni,EvtComplex c1,EvtComplex c2,EvtComplex c3,EvtComplex c4,EvtComplex c5,EvtComplex c6,double ialpha=137.036){
  EvtComplex i1(0,1);
  EvtComplex y(0,0);
  int i;
  for(i=0;i<ni;i++) y+=ki[i]*Gi[i]*Mi[i]/(q2-Mi[i]*Mi[i]-i1*Mi[i]*Gi[i]);
  return (-3.*ialpha*ialpha*y*EvtConst::pi*(3.*c1+c2+3.*c3+c4+3.*c5+c6));
}
