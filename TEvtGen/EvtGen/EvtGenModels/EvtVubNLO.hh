//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtVubNLO.hh
//
// Description:
// Class to generate inclusive B to X_u l nu decays according to various
// decay models. Implemtented are ACCM, parton-model and a QCD model.
//
// Modification history:
//
//   Sven Menke     January 17, 2001         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVUBNLO_HH
#define EVTVUBNLO_HH

#include <vector>
#include "EvtGenBase/EvtDecayIncoherent.hh"

class EvtParticle;
class RandGeneral;

class EvtVubNLO:public  EvtDecayIncoherent  {

public:
  
  EvtVubNLO() {}
  virtual ~EvtVubNLO();

  std::string getName();

  EvtDecayBase* clone();

  void initProbMax();

  void init();

  void decay(EvtParticle *p); 


private:

  // cache
  double _lbar;
  double _mupi2;

  double _mb;     // the b-quark pole mass in GeV 
  double _mB;
  double _lambdaSF;
  double _b;      // Parameter for the Fermi Motion 
  double _kpar;  
  double _mui; // renormalization scale (preferred value=1.5 GeV)
  double _SFNorm; // SF normalization 
  double _dGMax;  // max dGamma*p2 value;
  int    _nbins;
  int    _idSF;// which shape function?
  double * _masses;
  double * _weights;

  double _gmax;
  int _ngood,_ntot;


  double tripleDiff(double pp, double pl, double pm);
  double SFNorm(const std::vector<double> &coeffs);
  static double integrand(double omega, const std::vector<double> &coeffs);
  double F10(const std::vector<double> &coeffs);
  static double F1Int(double omega,const std::vector<double> &coeffs);
  double F20(const std::vector<double> &coeffs);
  static double F2Int(double omega,const std::vector<double> &coeffs);
  double F30(const std::vector<double> &coeffs);
  static double F3Int(double omega,const std::vector<double> &coeffs);
  static double g1(double y, double z);
  static double g2(double y, double z);
  static double g3(double y, double z);

  static double Gamma(double z);// Euler Gamma Function
  static double dgamma(double t, const std::vector<double> &c){  return pow(t,c[0]-1)*exp(-t);}
  static double Gamma(double z, double tmax);

  // theory parameters
  inline double mu_i(){return _mui;} // intermediate scale
  inline double mu_bar(){return _mui;} 
  inline double mu_h(){return _mb/sqrt(2.0);} // high scale
  inline double lambda1(){return -_mupi2;}
  
  // expansion coefficients for RGE
  static double beta0(int nf=4){return 11.-2./3.*nf;}
  static double beta1(int nf=4){return 34.*3.-38./3.*nf;}
  static double beta2(int nf=4){return 1428.5-5033./18.*nf+325./54.*nf*nf;}
  static double gamma0(){return 16./3.;}
  static double gamma1(int nf=4){return 4./3.*(49.85498-40./9.*nf);}
  static double gamma2(int nf=4){return 64./3.*(55.07242-8.58691*nf-nf*nf/27.);} /*  zeta3=1.20206 */
  static double gammap0(){return -20./3.;}
  static double gammap1(int nf=4){return -32./3.*(6.92653-0.9899*nf);} /* ??  zeta3=1.202 */


  // running constants

  static double alphas(double mu) ; 
  static double C_F(double mu){return (4.0/3.0)*alphas(mu)/4./EvtConst::pi;}

  // Shape Functions

  inline double lambda_SF(){ return _lambdaSF;}
  double lambda_bar(double omega0);
  inline double lambda2(){return 0.12;}
  double mu_pi2(double omega0);
  inline double lambda(double){ return _mB-_mb;}

  // specail for gaussian SF
  static double cGaus(double b){return pow(Gamma(1+b/2.)/Gamma((1+b)/2.),2);}

  double M0(double mui,double omega0);
  static double shapeFunction(double omega, const std::vector<double> &coeffs);
  static double expShapeFunction(double omega, const std::vector<double> &coeffs);
  static double gausShapeFunction(double omega, const std::vector<double> &coeffs);
  // SSF (not yet implemented)
  double subS(const std::vector<double> &coeffs );
  double subT(const std::vector<double> &coeffs);
  double subU(const std::vector<double> &coeffs);
  double subV(const std::vector<double> &coeffs);


  // Sudakov

  inline double S0(double a, double r){return -gamma0()/4/a/pow(beta0(),2)*(1/r-1+log(r));}
  inline double S1(double /*a*/, double r){return gamma0()/4./pow(beta0(),2)*(
								  pow(log(r),2)*beta1()/2./beta0()+(gamma1()/gamma0()-beta1()/beta0())*(1.-r+log(r))
								  );}
  inline double S2(double a, double r){return gamma0()*a/4./pow(beta0(),2)*(
									   -0.5*pow((1-r),2)*(
											 pow(beta1()/beta0(),2)-beta2()/beta0()-beta1()/beta0()*gamma1()/gamma0()+gamma2()/gamma0()
											 )
									   +(pow(beta1()/beta0(),2)-beta2()/beta0())*(1-r)*log(r)
									   +(beta1()/beta0()*gamma1()/gamma0()-beta2()/beta0())*(1-r+r*log(r))
									   );}
  inline double dSudakovdepsi(double mu1, double mu2){return S2(alphas(mu1)/(4*EvtConst::pi),alphas(mu2)/alphas(mu1));}
  inline double Sudakov(double mu1, double mu2, double epsi=0){double fp(4*EvtConst::pi);return S0(alphas(mu1)/fp,alphas(mu2)/alphas(mu1))+S1(alphas(mu1)/fp,alphas(mu2)/alphas(mu1))+epsi*dSudakovdepsi(mu1,mu2);}

  // RG 
  inline double dGdepsi(double mu1, double mu2){return 1./8./EvtConst::pi*(alphas(mu2)-alphas(mu1))*(gamma1()/beta0()-beta1()*gamma0()/pow(beta0(),2));}
  inline double aGamma(double mu1, double mu2, double epsi=0){return gamma0()/2/beta0()*log(alphas(mu2)/alphas(mu1))+epsi*dGdepsi( mu1, mu2);}
  inline double dgpdepsi(double mu1, double mu2){return 1./8./EvtConst::pi*(alphas(mu2)-alphas(mu1))*(gammap1()/beta0()-beta1()*gammap0()/pow(beta0(),2));}
  inline double agammap(double mu1, double mu2, double epsi=0){return gammap0()/2/beta0()*log(alphas(mu2)/alphas(mu1))+epsi*dgpdepsi( mu1, mu2);}
  inline double U1(double mu1, double mu2, double epsi=0){return exp(2*(Sudakov(mu1,mu2,epsi)-agammap(mu1,mu2,epsi)-aGamma(mu1,mu2,epsi)*log(_mb/mu1)));}
  inline double U1lo(double mu1, double mu2){return U1(mu1,mu2);}
  inline double U1nlo(double mu1, double mu2){return U1(mu1,mu2)*(1+2*(dSudakovdepsi(mu1,mu2)-dgpdepsi( mu1, mu2)-log(_mb/mu1)*dGdepsi( mu1, mu2)));}
  inline double alo(double mu1, double mu2){return -2*aGamma(mu1,mu2);}
  inline double anlo(double mu1, double mu2){return -2*dGdepsi(mu1,mu2);}

};

#endif

