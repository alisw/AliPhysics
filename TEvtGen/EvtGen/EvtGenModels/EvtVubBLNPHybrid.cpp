
//////////////////////////////////////////////////////////////////////
//
// Module: EvtVubBLNPHybrid.cc
//
// Description: Modeled on Riccardo Faccini's EvtVubNLO module
//
// tripleDiff from BLNP's notebook (based on BLNP4, hep-ph/0504071)
//
//////////////////////////////////////////////////////////////////

#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtVubBLNPHybrid.hh"
#include <string>
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenModels/EvtItgSimpsonIntegrator.hh"
#include "EvtGenModels/EvtItgPtrFunction.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenModels/EvtPFermi.hh"

// For incomplete gamma function
#include "math.h"
#include "signal.h"
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

using std::cout;
using std::endl;

EvtVubBLNPHybrid::EvtVubBLNPHybrid() 
  : _noHybrid(false), _storeWhat(true),
    _nbins_mX(0), _nbins_q2(0), _nbins_El(0), _nbins(0),
    _masscut(0.28), _bins_mX(0), _bins_q2(0), _bins_El(0),
    _weights(0)
{}


EvtVubBLNPHybrid::~EvtVubBLNPHybrid() {
  delete [] _bins_mX;
  delete [] _bins_q2;
  delete [] _bins_El;
  delete [] _weights;
}

std::string EvtVubBLNPHybrid::getName(){
  return "VUB_BLNPHYBRID";
}

EvtDecayBase *EvtVubBLNPHybrid::clone() {

  return new EvtVubBLNPHybrid;

}

void EvtVubBLNPHybrid::init() {
  
  // check that there are at least 3 arguments
  if (getNArg() < EvtVubBLNPHybrid::nParameters) {
    report(Severity::Error,"EvtVubBLNPHybrid") << "EvtVubBLNPHybrid generator expected "
				     << "at least " << EvtVubBLNPHybrid::nParameters
				     << " arguments but found: " << getNArg()
				     << "\nWill terminate execution!"<<endl;
    ::abort(); 
  } else if (getNArg() == EvtVubBLNPHybrid::nParameters) {
    report(Severity::Warning,"EvtVubBLNPHybrid") << "EvtVubBLNPHybrid: generate B -> Xu l nu events " 
				   << "without using the hybrid reweighting." 
				   << endl;
    _noHybrid = true;
  } else if (getNArg() < EvtVubBLNPHybrid::nParameters+EvtVubBLNPHybrid::nVariables) {
    report(Severity::Error,"EvtVubBLNPHybrid") << "EvtVubBLNPHybrid could not read number of bins for "
				     << "all variables used in the reweighting\n"
				     << "Will terminate execution!" << endl;
    ::abort();    
  }
  

  
  // get parameters (declared in the header file)
  
  // Input parameters
  mBB = 5.2792;
  lambda2 = 0.12;

  // Shape function parameters
  b = getArg(0);
  Lambda = getArg(1);
  Ecut = 1.8;
  wzero = mBB - 2*Ecut;

  // SF and SSF modes
  itype = (int)getArg(5);
  dtype = getArg(5);
  isubl = (int)getArg(6);

  // flags
  flag1 = (int)getArg(7);
  flag2 = (int)getArg(8);
  flag3 = (int)getArg(9);

  // Quark mass
  mb = 4.61;


  // hidden parameter what and SF stuff   
  const double xlow = 0;
  const double xhigh = mBB;
  const int aSize = 10000;
  EvtPFermi pFermi(Lambda,b);
  // pf is the cumulative distribution normalized to 1.
  _pf.resize(aSize);
  for(int i=0;i<aSize;i++){
    double what = xlow + (double)(i+0.5)/((double)aSize)*(xhigh-xlow);
    if ( i== 0 )
      _pf[i] = pFermi.getSFBLNP(what);
    else
      _pf[i] = _pf[i-1] + pFermi.getSFBLNP(what);
  }
  for (size_t i=0; i<_pf.size(); i++) {
    _pf[i]/=_pf[_pf.size()-1];
  }  



  // Matching scales
  muh = mBB*getArg(2); // 0.5
  mui = getArg(3); // 1.5
  mubar = getArg(4); // 1.5

  // Perturbative quantities
  CF = 4.0/3.0;
  CA = 3.0;
  double nf = 4.0;

  beta0 = 11.0/3.0*CA - 2.0/3.0*nf;
  beta1 = 34.0/3.0*CA*CA - 10.0/3.0*CA*nf - 2.0*CF*nf;
  beta2 = 2857.0/54.0*CA*CA*CA + (CF*CF - 205.0/18.0*CF*CA - 1415.0/54.0*CA*CA)*nf + (11.0/9.0*CF + 79.0/54.0*CA)*nf*nf;

  zeta3 = 1.0 + 1/8.0 + 1/27.0 + 1/64.0;

  Gamma0 = 4*CF;
  Gamma1 = CF*( (268.0/9.0 - 4.0*M_PI*M_PI/3.0)*CA - 40.0/9.0*nf);
  Gamma2 = 16*CF*( (245.0/24.0 - 67.0/54.0*M_PI*M_PI + + 11.0/180.0*pow(M_PI,4) + 11.0/6.0*zeta3)*CA*CA* + (-209.0/108.0 + 5.0/27.0*M_PI*M_PI - 7.0/3.0*zeta3)*CA*nf + (-55.0/24.0 + 2*zeta3)*CF*nf - nf*nf/27.0);

  gp0 = -5.0*CF;
  gp1 = -8.0*CF*( (3.0/16.0 - M_PI*M_PI/4.0 + 3*zeta3)*CF + (1549.0/432.0 + 7.0/48.0*M_PI*M_PI - 11.0/4.0*zeta3)*CA - (125.0/216.0 + M_PI*M_PI/24.0)*nf );

  // Lbar and mupisq

  Lbar = Lambda;  // all models
  mupisq = 3*Lambda*Lambda/b;
  if (itype == 1) mupisq = 3*Lambda*Lambda/b;
  if (itype == 2) mupisq = 3*Lambda*Lambda*(Gamma(1+0.5*b)*Gamma(0.5*b)/pow( Gamma(0.5 + 0.5*b), 2) - 1);

  // moment2 for SSFs
  moment2 = pow(0.3,3);

  // inputs for total rate (T for Total); use BLNP notebook defaults
  flagpower = 1;
  flag2loop = 1;

  // stuff for the integrator
  maxLoop = 20;
  //precision = 1.0e-3;
  precision = 2.0e-2;

  // vector of global variables, to pass to static functions (which can't access globals);
  gvars.push_back(0.0); // 0
  gvars.push_back(0.0); // 1
  gvars.push_back(mui); // 2
  gvars.push_back(b); // 3
  gvars.push_back(Lambda); // 4
  gvars.push_back(mBB); // 5
  gvars.push_back(mb); // 6
  gvars.push_back(wzero); // 7
  gvars.push_back(beta0); // 8
  gvars.push_back(beta1); // 9
  gvars.push_back(beta2); // 10
  gvars.push_back(dtype); // 11

  // check that there are 3 daughters and 10 arguments
  checkNDaug(3);
  // A. Volk: check for number of arguments is not necessary
  //checkNArg(10);
  
  if (_noHybrid) return;	// Without hybrid weighting, nothing else to do
 
  _nbins_mX = abs((int)getArg(10));
  _nbins_q2 = abs((int)getArg(11));
  _nbins_El = abs((int)getArg(12));

  int nextArg = EvtVubBLNPHybrid::nParameters + EvtVubBLNPHybrid::nVariables;

  _nbins = _nbins_mX*_nbins_q2*_nbins_El;	// Binning of weight table

  int expectArgs = nextArg + _nbins_mX +_nbins_q2 + _nbins_El + _nbins;
  
  if (getNArg() < expectArgs) {
    report(Severity::Error,"EvtVubBLNPHybrid")
      << " finds " << getNArg() << " arguments, expected " << expectArgs
      << ".  Something is wrong with the tables of weights or thresholds."
      << "\nWill terminate execution!" << endl;
    ::abort();        
  }

  // read bin boundaries from decay.dec
  int i;

  _bins_mX  = new double[_nbins_mX];
  for (i = 0; i < _nbins_mX; i++,nextArg++) {
    _bins_mX[i] = getArg(nextArg);
  }
  _masscut = _bins_mX[0];
  
  _bins_q2  = new double[_nbins_q2];
  for (i = 0; i < _nbins_q2; i++,nextArg++) {
    _bins_q2[i] = getArg(nextArg);    
  }
  
  _bins_El  = new double[_nbins_El];
  for (i = 0; i < _nbins_El; i++,nextArg++) {
    _bins_El[i] = getArg(nextArg);    
  }
  
  // read in weights (and rescale to range 0..1)
  readWeights(nextArg); 
  
}

void EvtVubBLNPHybrid::initProbMax() {
  noProbMax();
}

void EvtVubBLNPHybrid::decay(EvtParticle *Bmeson) {

  int j;
  
  EvtParticle *xuhad, *lepton, *neutrino;
  EvtVector4R p4;
  double Pp, Pm, Pl, pdf, EX, sh, qsq, El, ml, mpi, ratemax;
  
  double xhigh, xlow, what;
  double mX;
 
  
  bool rew(true);
  while(rew){

    Bmeson->initializePhaseSpace(getNDaug(), getDaugs());
    
    xuhad = Bmeson->getDaug(0);
    lepton = Bmeson->getDaug(1);
    neutrino = Bmeson ->getDaug(2);

    mBB = Bmeson->mass();
    ml = lepton->mass();

  
  
    //  get SF value 
    xlow = 0;
    xhigh = mBB;    
    // the case for alphas = 0 is not considered 
    what = 2*xhigh;
    while( what > xhigh || what < xlow ) {
      what = findBLNPWhat(); 
      what = xlow + what*(xhigh-xlow);
    }
  
  
  
    bool tryit = true;
  
    while (tryit) {
      
      // generate pp between 0 and 
      // Flat(min, max) gives R(max - min) + min, where R = random btwn 0 and 1

      Pp = EvtRandom::Flat(0, mBB); // P+ = EX - |PX|
      Pl = EvtRandom::Flat(0, mBB);  // mBB - 2El
      Pm = EvtRandom::Flat(0, mBB); // P- = EX + |PX|

      sh = Pm*Pp;
      EX = 0.5*(Pm + Pp);
      qsq = (mBB - Pp)*(mBB - Pm);
      El = 0.5*(mBB - Pl);
      
      // Need maximum rate.  Waiting for Mr. Paz to give it to me. 
      // Meanwhile, use this.
      ratemax = 3.0;  // From trial and error - most events below 3.0

      // kinematic bounds (Eq. 2)
      mpi = 0.14;
      if ((Pp > 0)&&(Pp <= Pl)&&(Pl <= Pm)&&(Pm < mBB)&&(El > ml)&&(sh > 4*mpi*mpi)) {
	
	// Probability of pass proportional to PDF
	pdf = rate3(Pp, Pl, Pm);
	double testRan = EvtRandom::Flat(0., ratemax);
	if (pdf >= testRan) tryit = false;
      }
    }
    
    // compute all kinematic variables needed for reweighting
    mX = sqrt(sh);
    
    // Reweighting in bins of mX, q2, El 
    if (_nbins>0) {
      double xran1 = EvtRandom::Flat();
      double w = 1.0;
      if (!_noHybrid) w = getWeight(mX, qsq, El); 
      if ( w >= xran1 ) rew = false;
    } 
    else {
      rew = false;
    }
  }
  // o.k. we have the three kineamtic variables 
  // now calculate a flat cos Theta_H [-1,1] distribution of the 
  // hadron flight direction w.r.t the B flight direction 
  // because the B is a scalar and should decay isotropic.
  // Then chose a flat Phi_H [0,2Pi] w.r.t the B flight direction 
  // and and a flat Phi_L [0,2Pi] in the W restframe w.r.t the 
  // W flight direction.
  
  double ctH = EvtRandom::Flat(-1,1);
  double phH = EvtRandom::Flat(0,2*M_PI);
  double phL = EvtRandom::Flat(0,2*M_PI);

  // now compute the four vectors in the B Meson restframe
    
  double ptmp,sttmp;
  // calculate the hadron 4 vector in the B Meson restframe
  
  sttmp = sqrt(1-ctH*ctH);
  ptmp = sqrt(EX*EX-sh);
  double pHB[4] = {EX,ptmp*sttmp*cos(phH),ptmp*sttmp*sin(phH),ptmp*ctH};
  p4.set(pHB[0],pHB[1],pHB[2],pHB[3]);
  xuhad->init( getDaug(0), p4);
  

  if (_storeWhat ) {
    // cludge to store the hidden parameter what with the decay; 
    // the lifetime of the Xu is abused for this purpose.
    // tau = 1 ps corresponds to ctau = 0.3 mm -> in order to
    // stay well below BaBars sensitivity we take what/(10000 GeV).
    // To extract what back from the StdHepTrk its necessary to get
    // delta_ctau = Xu->decayVtx()->point().distanceTo(XuDaughter->decayVtx()->point());
    //
    // what = delta_ctau * 100000 * Mass_Xu/Momentum_Xu     
    //
    xuhad->setLifetime(what/10000.);
  }
  
  
  // calculate the W 4 vector in the B Meson restrframe

  double apWB = ptmp;
  double pWB[4] = {mBB-EX,-pHB[1],-pHB[2],-pHB[3]};

  // first go in the W restframe and calculate the lepton and
  // the neutrino in the W frame

  double mW2   = mBB*mBB + sh - 2*mBB*EX;
  double beta  = ptmp/pWB[0];
  double gamma = pWB[0]/sqrt(mW2);

  double pLW[4];
    
  ptmp = (mW2-ml*ml)/2/sqrt(mW2);
  pLW[0] = sqrt(ml*ml + ptmp*ptmp);

  double ctL = (El - gamma*pLW[0])/beta/gamma/ptmp;
  if ( ctL < -1 ) ctL = -1;
  if ( ctL >  1 ) ctL =  1;
  sttmp = sqrt(1-ctL*ctL);

  // eX' = eZ x eW
  double xW[3] = {-pWB[2],pWB[1],0};
  // eZ' = eW
  double zW[3] = {pWB[1]/apWB,pWB[2]/apWB,pWB[3]/apWB};
  
  double lx = sqrt(xW[0]*xW[0]+xW[1]*xW[1]);
  for (j=0;j<2;j++) 
    xW[j] /= lx;

  // eY' = eZ' x eX'
  double yW[3] = {-pWB[1]*pWB[3],-pWB[2]*pWB[3],pWB[1]*pWB[1]+pWB[2]*pWB[2]};
  double ly = sqrt(yW[0]*yW[0]+yW[1]*yW[1]+yW[2]*yW[2]);
  for (j=0;j<3;j++) 
    yW[j] /= ly;

  // p_lep = |p_lep| * (  sin(Theta) * cos(Phi) * eX'
  //                    + sin(Theta) * sin(Phi) * eY'
  //                    + cos(Theta) *            eZ')
  for (j=0;j<3;j++)
    pLW[j+1] = sttmp*cos(phL)*ptmp*xW[j] 
      +        sttmp*sin(phL)*ptmp*yW[j]
      +          ctL         *ptmp*zW[j];

  double apLW = ptmp;
    
  // boost them back in the B Meson restframe
  
  double appLB = beta*gamma*pLW[0] + gamma*ctL*apLW;
 
  ptmp = sqrt(El*El-ml*ml);
  double ctLL = appLB/ptmp;

  if ( ctLL >  1 ) ctLL =  1;
  if ( ctLL < -1 ) ctLL = -1;
    
  double pLB[4] = {El,0,0,0};
  double pNB[4] = {pWB[0]-El,0,0,0};

  for (j=1;j<4;j++) {
    pLB[j] = pLW[j] + (ctLL*ptmp - ctL*apLW)/apWB*pWB[j];
    pNB[j] = pWB[j] - pLB[j];
  }

  p4.set(pLB[0],pLB[1],pLB[2],pLB[3]);
  lepton->init( getDaug(1), p4);

  p4.set(pNB[0],pNB[1],pNB[2],pNB[3]);
  neutrino->init( getDaug(2), p4);

  return ;

}

double EvtVubBLNPHybrid::rate3(double Pp, double Pl, double Pm) {

  // rate3 in units of GF^2*Vub^2/pi^3

  double factor = 1.0/16*(mBB-Pp)*U1lo(muh, mui)*pow( (Pm - Pp)/(mBB - Pp), alo(muh, mui));

  double doneJS = DoneJS(Pp, Pm, mui);
  double done1 = Done1(Pp, Pm, mui);
  double done2 = Done2(Pp, Pm, mui);
  double done3 = Done3(Pp, Pm, mui);

  // The EvtSimpsonIntegrator returns zero for bad integrals.
  // So if any of the integrals are zero (ie bad), return zero.
  // This will cause pdf = 0, so the event will not pass.
  // I hope this will not introduce a bias.
  if (doneJS*done1*done2*done3 == 0.0) {
    //cout << "Integral failed: (Pp, Pm, Pl) = (" << Pp << ", " << Pm << ", " << Pl << ")" << endl;
    return 0.0;
  }
  //  if (doneJS*done1*done2*done3 != 0.0) {
  //    cout << "Integral OK: (Pp, Pm, Pl) = (" << Pp << ", " << Pm << ", " << Pl << ")" << endl;
  //}

  double f1 = F1(Pp, Pm, muh, mui, mubar, doneJS, done1);
  double f2 = F2(Pp, Pm, muh, mui, mubar, done3);
  double f3 = F3(Pp, Pm, muh, mui, mubar, done2);
  double answer = factor*( (mBB + Pl - Pp - Pm)*(Pm - Pl)*f1 + 2*(Pl - Pp)*(Pm - Pl)*f2 + (mBB - Pm)*(Pm - Pp)*f3 );
  return answer;

}

double EvtVubBLNPHybrid::F1(double Pp, double Pm, double muh, double mui, double mubar, double doneJS, double done1) {

  std::vector<double> vars(12);
  vars[0] = Pp;
  vars[1] = Pm;
  for (int j=2;j<12;j++) {vars[j] = gvars[j];}

  double y = (Pm - Pp)/(mBB - Pp);
  double ah = CF*alphas(muh, vars)/4/M_PI;
  double ai = CF*alphas(mui, vars)/4/M_PI;
  double abar = CF*alphas(mubar, vars)/4/M_PI;
  double lambda1 = -mupisq;

  double t1 = -4*ai/(Pp - Lbar)*(2*log((Pp - Lbar)/mui) + 1);
  double t2 = 1 + dU1nlo(muh, mui) + anlo(muh, mui)*log(y);
  double t3 = -4.0*pow(log(y*mb/muh),2) + 10.0*log(y*mb/muh) - 4.0*log(y) - 2.0*log(y)/(1-y) - 4.0*PolyLog(2, 1-y) - M_PI*M_PI/6.0 - 12.0;
  double t4 = 2*pow( log(y*mb*Pp/(mui*mui)), 2) - 3*log(y*mb*Pp/(mui*mui)) + 7 - M_PI*M_PI;

  double t5 = -wS(Pp) + 2*t(Pp) + (1.0/y - 1.0)*(u(Pp) - v(Pp));
  double t6 = -(lambda1 + 3.0*lambda2)/3.0 + 1.0/pow(y,2)*(4.0/3.0*lambda1 - 2.0*lambda2);

  double shapePp = Shat(Pp, vars);

  double answer = (t2 + ah*t3 + ai*t4)*shapePp + ai*doneJS + 1/(mBB - Pp)*(flag2*abar*done1 + flag1*t5) + 1/pow(mBB - Pp, 2)*flag3*shapePp*t6;
  if (Pp > Lbar + mui/exp(0.5)) answer = answer + t1;
  return answer;

}

double EvtVubBLNPHybrid::F2(double Pp, double Pm, double muh, double /* mui */, double mubar, double done3) {
  
  std::vector<double> vars(12);
  vars[0] = Pp;
  vars[1] = Pm;
  for (int j=2;j<12;j++) {vars[j] = gvars[j];}

  double y = (Pm - Pp)/(mBB - Pp);
  double lambda1 = -mupisq;
  double ah = CF*alphas(muh, vars)/4/M_PI;
  double abar = CF*alphas(mubar, vars)/4/M_PI;

  double t6 = -wS(Pp) - 2*t(Pp) + 1.0/y*(t(Pp) + v(Pp));
  double t7 = 1/pow(y,2)*(2.0/3.0*lambda1 + 4.0*lambda2) - 1/y*(2.0/3.0*lambda1 + 3.0/2.0*lambda2);

  double shapePp = Shat(Pp, vars);

  double answer = ah*log(y)/(1-y)*shapePp + 1/(mBB - Pp)*(flag2*abar*0.5*done3 + flag1/y*t6) + 1.0/pow(mBB - Pp,2)*flag3*shapePp*t7;
  return answer;

}

double EvtVubBLNPHybrid::F3(double Pp, double Pm, double /*muh*/, double /* mui */, double mubar, double done2) {

  std::vector<double> vars(12);
  vars[0] = Pp;
  vars[1] = Pm;
  for (int j=2;j<12;j++) {vars[j] = gvars[j];}
  
  double y = (Pm - Pp)/(mBB - Pp);
  double lambda1 = -mupisq;
  double abar = CF*alphas(mubar, vars)/4/M_PI;

  double t7 = 1.0/pow(y,2)*(-2.0/3.0*lambda1 + lambda2);

  double shapePp = Shat(Pp, vars);

  double answer = 1.0/(Pm - Pp)*flag2*0.5*y*abar*done2 + 1.0/pow(mBB-Pp,2)*flag3*shapePp*t7;
  return answer;

}

double EvtVubBLNPHybrid::DoneJS(double Pp, double Pm, double /* mui */) {

  std::vector<double> vars(12);
  vars[0] = Pp;
  vars[1] = Pm;
  for (int j=2;j<12;j++) {vars[j] = gvars[j];}
  
  double lowerlim = 0.001*Pp;
  double upperlim = (1.0-0.001)*Pp;

  EvtItgPtrFunction *func = new EvtItgPtrFunction(&IntJS, lowerlim, upperlim, vars);
  EvtItgSimpsonIntegrator *integ = new EvtItgSimpsonIntegrator(*func, precision, maxLoop);
  double myintegral = integ->evaluate(lowerlim, upperlim);
  delete integ;
  delete func;
  return myintegral;

}

double EvtVubBLNPHybrid::Done1(double Pp, double Pm, double /* mui */) {

  std::vector<double> vars(12);
  vars[0] = Pp;
  vars[1] = Pm;
  for (int j=2;j<12;j++) {vars[j] = gvars[j];}

  double lowerlim = 0.001*Pp;
  double upperlim = (1.0-0.001)*Pp;

  EvtItgPtrFunction *func = new EvtItgPtrFunction(&Int1, lowerlim, upperlim, vars);
  EvtItgSimpsonIntegrator *integ = new EvtItgSimpsonIntegrator(*func, precision, maxLoop);
  double myintegral = integ->evaluate(lowerlim, upperlim);
  delete integ;
  delete func;
  return myintegral;

}

double EvtVubBLNPHybrid::Done2(double Pp, double Pm, double /* mui */ ) {

  std::vector<double> vars(12);
  vars[0] = Pp;
  vars[1] = Pm;
  for (int j=2;j<12;j++) {vars[j] = gvars[j];}

  double lowerlim = 0.001*Pp;
  double upperlim = (1.0-0.001)*Pp;

  EvtItgPtrFunction *func = new EvtItgPtrFunction(&Int2, lowerlim, upperlim, vars);
  EvtItgSimpsonIntegrator *integ = new EvtItgSimpsonIntegrator(*func, precision, maxLoop);
  double myintegral = integ->evaluate(lowerlim, upperlim);
  delete integ;
  delete func;
  return myintegral;

}

double EvtVubBLNPHybrid::Done3(double Pp, double Pm, double /* mui */) {

  std::vector<double> vars(12);
  vars[0] = Pp;
  vars[1] = Pm;
  for (int j=2;j<12;j++) {vars[j] = gvars[j];}

  double lowerlim = 0.001*Pp;
  double upperlim = (1.0-0.001)*Pp;  

  EvtItgPtrFunction *func = new EvtItgPtrFunction(&Int3, lowerlim, upperlim, vars);
  EvtItgSimpsonIntegrator *integ = new EvtItgSimpsonIntegrator(*func, precision, maxLoop);
  double myintegral = integ->evaluate(lowerlim, upperlim);
  delete integ;
  delete func;
  return myintegral;

}

double EvtVubBLNPHybrid::Int1(double what, const std::vector<double> &vars) {
  return Shat(what, vars)*g1(what, vars);
}

double EvtVubBLNPHybrid::Int2(double what, const std::vector<double> &vars) {
  return Shat(what, vars)*g2(what, vars);
}

double EvtVubBLNPHybrid::Int3(double what, const std::vector<double> &vars) {
  return Shat(what, vars)*g3(what, vars);
}

double EvtVubBLNPHybrid::IntJS(double what, const std::vector<double> &vars) {
  
  double Pp = vars[0];
  double Pm = vars[1];
  double mui = vars[2];
  double mBB = vars[5];
  double mb = vars[6];
  double y = (Pm - Pp)/(mBB - Pp);
  
  return 1/(Pp-what)*(Shat(what, vars) - Shat(Pp, vars))*(4*log(y*mb*(Pp-what)/(mui*mui)) - 3);
}

double EvtVubBLNPHybrid::g1(double w, const std::vector<double> &vars) {

  double Pp = vars[0];
  double Pm = vars[1];
  double mBB = vars[5];
  double y = (Pm - Pp)/(mBB - Pp);
  double x = (Pp - w)/(mBB - Pp);

  double q1 = (1+x)*(1+x)*y*(x+y);
  double q2 = y*(-9 + 10*y) + x*x*(-12.0 + 13.0*y) + 2*x*(-8.0 + 6*y + 3*y*y);
  double q3 = 4/x*log(y + y/x);
  double q4 = 3.0*pow(x,4)*(-2.0 + y) - 2*pow(y,3) - 4*pow(x,3)*(2.0+y) - 2*x*y*y*(4+y) - x*x*y*(12 + 4*y + y*y);
  double q5 = log(1 + y/x);

  double answer = q2/q1 - q3 - 2*q4*q5/(q1*y*x);
  return answer;

}

double EvtVubBLNPHybrid::g2(double w, const std::vector<double> &vars) {

  double Pp = vars[0];
  double Pm = vars[1];
  double mBB = vars[5];
  double y = (Pm - Pp)/(mBB - Pp);
  double x = (Pp - w)/(mBB - Pp);

  double q1 = (1+x)*(1+x)*pow(y,3)*(x+y);
  double q2 = 10.0*pow(x,4) + y*y + 3.0*pow(x,2)*y*(10.0+y) + pow(x,3)*(12.0+19.0*y) + x*y*(8.0 + 4.0*y + y*y);
  double q3 = 5*pow(x,4) + 2.0*y*y + 6.0*pow(x,3)*(1.0+2.0*y) + 4.0*x*y*(1+2.0*y) + x*x*y*(18.0+5.0*y);
  double q4 = log(1 + y/x);

  double answer = 2.0/q1*( y*q2 - 2*x*q3*q4);
  return answer;

}

double EvtVubBLNPHybrid::g3(double w, const std::vector<double> &vars) {

  double Pp = vars[0];
  double Pm = vars[1];
  double mBB = vars[5];
  double y = (Pm - Pp)/(mBB - Pp);
  double x = (Pp - w)/(mBB - Pp);

  double q1 = (1+x)*(1+x)*pow(y,3)*(x+y);
  double q2 =  2.0*pow(y,3)*(-11.0+2.0*y) - 10.0*pow(x,4)*(6 - 6*y + y*y) + x*y*y*(-94.0 + 29.0*y + 2.0*y*y) + 2.0*x*x*y*(-72.0 +18.0*y + 13.0*y*y) - x*x*x*(72.0 + 42.0*y - 70.0*y*y + 3.0*y*y*y);
  double q3 =  -6.0*x*(-5.0+y)*pow(y,3) + 4*pow(y,4) + 5*pow(x,5)*(6-6*y + y*y) - 4*x*x*y*y*(-20.0 + 6*y + y*y) + pow(x,3)*y*(90.0 - 10.0*y - 28.0*y*y + y*y*y) + pow(x,4)*(36.0 + 36.0*y - 50.0*y*y + 4*y*y*y);
  double q4 = log(1 + y/x);

  double answer = q2/q1 + 2/q1/y*q3*q4;
  return answer;

}


double EvtVubBLNPHybrid::Shat(double w, const std::vector<double> &vars) {

  double mui = vars[2];
  double b = vars[3];
  double Lambda = vars[4];
  double wzero = vars[7];
  int itype = (int)vars[11];

  double norm = 0.0;
  double shape = 0.0;

  if (itype == 1) {

    double Lambar = (Lambda/b)*(Gamma(1+b)-Gamma(1+b,b*wzero/Lambda))/(Gamma(b) - Gamma(b, b*wzero/Lambda));
    double muf = wzero - Lambar;
    double mupisq = 3*pow(Lambda,2)/pow(b,2)*(Gamma(2+b) - Gamma(2+b, b*wzero/Lambda))/(Gamma(b) - Gamma(b, b*wzero/Lambda)) - 3*Lambar*Lambar;
    norm = Mzero(muf, mui, mupisq, vars)*Gamma(b)/(Gamma(b) - Gamma(b, b*wzero/Lambda));
    shape = pow(b,b)/Lambda/Gamma(b)*pow(w/Lambda, b-1)*exp(-b*w/Lambda);
  }

  if (itype == 2) {
    double dcoef = pow( Gamma(0.5*(1+b))/Gamma(0.5*b), 2);
    double t1 =  wzero*wzero*dcoef/(Lambda*Lambda);
    double Lambar = Lambda*(Gamma(0.5*(1+b)) - Gamma(0.5*(1+b),t1))/pow(dcoef, 0.5)/(Gamma(0.5*b) - Gamma(0.5*b, t1));
    double muf = wzero - Lambar;
    double mupisq = 3*Lambda*Lambda*( Gamma(1+0.5*b) - Gamma(1+0.5*b, t1))/dcoef/(Gamma(0.5*b) - Gamma(0.5*b, t1)) - 3*Lambar*Lambar;
    norm = Mzero(muf, mui, mupisq, vars)*Gamma(0.5*b)/(Gamma(0.5*b) - Gamma(0.5*b, wzero*wzero*dcoef/(Lambda*Lambda)));
    shape = 2*pow(dcoef, 0.5*b)/Lambda/Gamma(0.5*b)*pow(w/Lambda, b-1)*exp(-dcoef*w*w/(Lambda*Lambda));
  }

  double answer = norm*shape;
  return answer;
}

double EvtVubBLNPHybrid::Mzero(double muf, double mu, double mupisq, const std::vector<double> &vars) {

  double CF = 4.0/3.0;
  double amu = CF*alphas(mu, vars)/M_PI;
  double answer = 1 - amu*( pow(log(muf/mu), 2) + log(muf/mu) + M_PI*M_PI/24.0) + amu*(log(muf/mu) - 0.5)*mupisq/(3*muf*muf);
  return answer;

}

double EvtVubBLNPHybrid::wS(double w) {

  double answer = (Lbar - w)*Shat(w, gvars);
  return answer;
}

double EvtVubBLNPHybrid::t(double w) {

  double t1 = -3*lambda2/mupisq*(Lbar - w)*Shat(w, gvars);
  double myf = myfunction(w, Lbar, moment2);
  double myBIK = myfunctionBIK(w, Lbar, moment2);
  double answer = t1;

  if (isubl == 1) answer = t1;
  if (isubl == 3) answer = t1 - myf;
  if (isubl == 4) answer = t1 + myf;
  if (isubl == 5) answer = t1 - myBIK;
  if (isubl == 6) answer = t1 + myBIK;

  return answer;
}

double EvtVubBLNPHybrid::u(double w) {

  double u1 = -2*(Lbar - w)*Shat(w, gvars);
  double myf = myfunction(w, Lbar, moment2);
  double myBIK = myfunctionBIK(w, Lbar, moment2);
  double answer = u1;

  if (isubl == 1) answer = u1;
  if (isubl == 3) answer = u1 + myf;
  if (isubl == 4) answer = u1 - myf;
  if (isubl == 5) answer = u1 + myBIK;
  if (isubl == 6) answer = u1 - myBIK;

  return answer;
}

double EvtVubBLNPHybrid::v(double w) {

  double v1 = 3*lambda2/mupisq*(Lbar - w)*Shat(w, gvars);
  double myf = myfunction(w, Lbar, moment2);
  double myBIK = myfunctionBIK(w, Lbar, moment2);
  double answer = v1;

  if (isubl == 1) answer = v1;
  if (isubl == 3) answer = v1 - myf;
  if (isubl == 4) answer = v1 + myf;
  if (isubl == 5) answer = v1 - myBIK;
  if (isubl == 6) answer = v1 + myBIK;

  return answer;
}

double EvtVubBLNPHybrid::myfunction(double w, double Lbar, double mom2) {

  double bval = 5.0;
  double x = w/Lbar;
  double factor = 0.5*mom2*pow(bval/Lbar, 3);
  double answer = factor*exp(-bval*x)*(1 - 2*bval*x + 0.5*bval*bval*x*x);
  return answer;

}

double EvtVubBLNPHybrid::myfunctionBIK(double w, double Lbar, double /* mom2 */) {

  double aval = 10.0;
  double normBIK = (4 - M_PI)*M_PI*M_PI/8/(2-M_PI)/aval + 1;
  double z = 3*M_PI*w/8/Lbar;
  double q = M_PI*M_PI*2*pow(M_PI*aval, 0.5)*exp(-aval*z*z)/(4*M_PI - 8)*(1 - 2*pow(aval/M_PI, 0.5)*z) + 8/pow(1+z*z, 4)*(z*log(z) + 0.5*z*(1+z*z) - M_PI/4*(1-z*z));
  double answer = q/normBIK;
  return answer;

}

double EvtVubBLNPHybrid::dU1nlo(double muh, double mui) { 

  double ai = alphas(mui, gvars);
  double ah = alphas(muh, gvars);

  double q1 = (ah - ai)/(4*M_PI*beta0);
  double q2 = log(mb/muh)*Gamma1 + gp1;
  double q3 = 4*beta1*(log(mb/muh)*Gamma0 + gp0) + Gamma2*(1-ai/ah);
  double q4 = beta1*beta1*Gamma0*(-1.0 + ai/ah)/(4*pow(beta0,3));
  double q5 = -beta2*Gamma0*(1.0 + ai/ah) + beta1*Gamma1*(3 - ai/ah);
  double q6 = beta1*beta1*Gamma0*(ah - ai)/beta0 - beta2*Gamma0*ah + beta1*Gamma1*ai;
  
  double answer = q1*(q2 - q3/4/beta0 + q4 + q5/(4*beta0*beta0)) + 1/(8*M_PI*beta0*beta0*beta0)*log(ai/ah)*q6;
  return answer;
}

double EvtVubBLNPHybrid::U1lo(double muh, double mui) {
  double epsilon = 0.0;
  double answer = pow(mb/muh, -2*aGamma(muh, mui, epsilon))*exp(2*Sfun(muh, mui, epsilon) - 2*agp(muh, mui, epsilon));
  return answer;
}

double EvtVubBLNPHybrid::Sfun(double mu1, double mu2, double epsilon) {
  double a1 = alphas(mu1, gvars)/4/M_PI;
  double a2 = alphas(mu2, gvars)/alphas(mu1, gvars);

  double answer = S0(a1,a2) + S1(a1,a2) + epsilon*S2(a1,a2);
  return answer;

}

double EvtVubBLNPHybrid::S0(double a1, double r) {
  double answer = -Gamma0/(4.0*beta0*beta0*a1)*(-1.0 + 1.0/r + log(r));
  return answer;
}

double EvtVubBLNPHybrid::S1(double /* a1 */ , double r) {
  double answer = Gamma0/(4*beta0*beta0)*(0.5*log(r)*log(r)*beta1/beta0 + (Gamma1/Gamma0 - beta1/beta0)*(1 - r + log(r)));
  return answer;
}

double EvtVubBLNPHybrid::S2(double a1, double r) {

  double w1 = pow(beta1,2)/pow(beta0,2) - beta2/beta0 - beta1*Gamma1/(beta0*Gamma0) + Gamma2/Gamma0;
  double w2 = pow(beta1,2)/pow(beta0,2) - beta2/beta0;
  double w3 = beta1*Gamma1/(beta0*Gamma0) - beta2/beta0;
  double w4 = a1*Gamma0/(4*beta0*beta0);

  double answer = w4*(-0.5*pow(1-r,2)*w1 + w2*(1-r)*log(r) + w3*(1-r+r*log(r)));
  return answer;
}

double EvtVubBLNPHybrid::aGamma(double mu1, double mu2, double epsilon) {
  double a1 = alphas(mu1, gvars);
  double a2 = alphas(mu2, gvars);
  double answer = Gamma0/(2*beta0)*log(a2/a1) + epsilon*(a2-a1)/(8.0*M_PI)*(Gamma1/beta0 - beta1*Gamma0/(beta0*beta0));
  return answer;
}

double EvtVubBLNPHybrid::agp(double mu1, double mu2, double epsilon) { 
  double a1 = alphas(mu1, gvars);
  double a2 = alphas(mu2, gvars);
  double answer = gp0/(2*beta0)*log(a2/a1) + epsilon*(a2-a1)/(8.0*M_PI)*(gp1/beta0 - beta1*gp0/(beta0*beta0));
  return answer;
}

double EvtVubBLNPHybrid::alo(double muh, double mui) { return -2.0*aGamma(muh, mui, 0);}

double EvtVubBLNPHybrid::anlo(double muh, double mui) {   // d/depsilon of aGamma

  double ah = alphas(muh, gvars);
  double ai = alphas(mui, gvars);
  double answer = (ah-ai)/(8.0*M_PI)*(Gamma1/beta0 - beta1*Gamma0/(beta0*beta0));
  return answer;
}

double EvtVubBLNPHybrid::alphas(double mu, const std::vector<double> &vars) {

  // Note: Lambda4 and Lambda5 depend on mbMS = 4.25
  // So if you change mbMS, then you will have to recalculate them.

  double beta0 = vars[8];
  double beta1 = vars[9];
  double beta2 = vars[10];
  
  double Lambda4 = 0.298791;
  double lg = 2*log(mu/Lambda4);
  double answer = 4*M_PI/(beta0*lg)*( 1 - beta1*log(lg)/(beta0*beta0*lg) + beta1*beta1/(beta0*beta0*beta0*beta0*lg*lg)*( (log(lg) - 0.5)*(log(lg) - 0.5) - 5.0/4.0 + beta2*beta0/(beta1*beta1)));
  return answer;
    
}

double EvtVubBLNPHybrid::PolyLog(double v, double z) {

  if (z >= 1) cout << "Error in EvtVubBLNPHybrid: 2nd argument to PolyLog is >= 1." << endl;

  double sum = 0.0;
  for (int k=1; k<101; k++) { 
    sum = sum + pow(z,k)/pow(k,v);
  }
  return sum;
}

double EvtVubBLNPHybrid::Gamma(double z)
{
   if (z<=0) return 0;

   double v = lgamma(z);
   return exp(v);
}

double EvtVubBLNPHybrid::Gamma(double a, double x)
{
    double LogGamma;
    /*    if (x<0.0 || a<= 0.0) raise(SIGFPE);*/
    if(x<0.0) x=0.0;
    if(a<=0.0)a=1.e-50;
    LogGamma = lgamma(a);
    if (x < (a+1.0)) 
        return gamser(a,x,LogGamma);
    else 
        return 1.0-gammcf(a,x,LogGamma);
}

/* ------------------Incomplete gamma function-----------------*/
/* ------------------via its series representation-------------*/
              
double EvtVubBLNPHybrid::gamser(double a, double x, double LogGamma)
{
    double n;
    double ap,del,sum;

    ap=a;
    del=sum=1.0/a;
    for (n=1;n<ITMAX;n++) {
        ++ap;
        del *= x/ap;
        sum += del;
        if (fabs(del) < fabs(sum)*EPS) return sum*exp(-x + a*log(x) - LogGamma);
    }
    raise(SIGFPE);

    return 0.0;
}        

/* ------------------Incomplete gamma function complement------*/
/* ------------------via its continued fraction representation-*/

double EvtVubBLNPHybrid::gammcf(double a, double x, double LogGamma) {
  
    double an,b,c,d,del,h;
    int i;

    b = x + 1.0 -a;
    c = 1.0/FPMIN;
    d = 1.0/b;
    h = d;
    for (i=1;i<ITMAX;i++) {
        an = -i*(i-a);
        b+=2.0;
        d=an*d+b;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = b+an/c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS) return exp(-x+a*log(x)-LogGamma)*h;  
    }
    raise(SIGFPE);

    return 0.0;

}


double EvtVubBLNPHybrid::findBLNPWhat() {

  double ranNum=EvtRandom::Flat();
  double oOverBins= 1.0/(float(_pf.size()));
  int nBinsBelow = 0;	  // largest k such that I[k] is known to be <= rand
  int nBinsAbove = _pf.size();  // largest k such that I[k] is known to be >  rand
  int middle;
  
  while (nBinsAbove > nBinsBelow+1) {
    middle = (nBinsAbove + nBinsBelow+1)>>1;
    if (ranNum >= _pf[middle]) {
      nBinsBelow = middle;
    } else {
      nBinsAbove = middle;
    }
  } 

  double bSize = _pf[nBinsAbove] - _pf[nBinsBelow];
  // binMeasure is always aProbFunc[nBinsBelow], 
  
  if ( bSize == 0 ) { 
    // rand lies right in a bin of measure 0.  Simply return the center
    // of the range of that bin.  (Any value between k/N and (k+1)/N is 
    // equally good, in this rare case.)
    return (nBinsBelow + .5) * oOverBins;
  }
  
  double bFract = (ranNum - _pf[nBinsBelow]) / bSize;
  
  return (nBinsBelow + bFract) * oOverBins;
  
} 

double EvtVubBLNPHybrid::getWeight(double mX, double q2, double El) {

  int ibin_mX = -1;
  int ibin_q2 = -1;
  int ibin_El = -1;

  for (int i = 0; i < _nbins_mX; i++) {
    if (mX >= _bins_mX[i]) ibin_mX = i;
  }
  for (int i = 0; i < _nbins_q2; i++) {
    if (q2 >= _bins_q2[i]) ibin_q2 = i;
  }
  for (int i = 0; i < _nbins_El; i++) {
    if (El >= _bins_El[i]) ibin_El = i;
  }
  int ibin = ibin_mX + ibin_q2*_nbins_mX + ibin_El*_nbins_mX*_nbins_q2;

  if ( (ibin_mX < 0) || (ibin_q2 < 0) || (ibin_El < 0) ) {
    report(Severity::Error,"EvtVubHybrid") << "Cannot determine hybrid weight "
                                 << "for this event " 
				 << "-> assign weight = 0" << endl;
    return 0.0;
  }

  return _weights[ibin];
}


void EvtVubBLNPHybrid::readWeights(int startArg) {
  _weights  = new double[_nbins];

  double maxw = 0.0;
  for (int i = 0; i < _nbins; i++, startArg++) {
    _weights[i] = getArg(startArg);
    if (_weights[i] > maxw) maxw = _weights[i];
  }

  if (maxw == 0) {
    report(Severity::Error,"EvtVubBLNPHybrid") << "EvtVub generator expected at least one " 
				     << " weight > 0, but found none! " 
				     << "Will terminate execution!"<<endl;
    ::abort();
  }

  // rescale weights (to be in range 0..1)
  for (int i = 0; i < _nbins; i++) {
    _weights[i] /= maxw;
  }
}
