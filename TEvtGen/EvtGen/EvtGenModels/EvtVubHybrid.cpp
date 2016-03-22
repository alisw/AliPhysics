//---------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtVubHybrid.cc
//
// Description: Routine to decay a particle according to phase space. 
//
// Modification history:
//
//   Jochen Dingfelder February 1, 2005  Created Module as update of the
//                                       original module EvtVub by Sven Menke 
//---------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtVubHybrid.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenModels/EvtPFermi.hh"
#include "EvtGenModels/EvtVubdGamma.hh"
#include "EvtGenBase/EvtRandom.hh"

#include <string>
#include <iostream>
#include <fstream>
using std::ifstream;
using std::cout;
using std::endl;


// _noHybrid will be set TRUE if the DECAY.DEC file has no binning or weights
// _storeQplus should alwasy be TRUE: writes out Fermi motion parameter

EvtVubHybrid::EvtVubHybrid() 
  : _noHybrid(false), _storeQplus(true),
    _mb(4.62), _a(2.27), _alphas(0.22), _dGMax(3.),
    _nbins_mX(0), _nbins_q2(0), _nbins_El(0), _nbins(0),
    _masscut(0.28), _bins_mX(0), _bins_q2(0), _bins_El(0),
    _weights(0), _dGamma(0)
{}

EvtVubHybrid::~EvtVubHybrid() {
  delete _dGamma;
  delete [] _bins_mX;
  delete [] _bins_q2;
  delete [] _bins_El;
  delete [] _weights;

}

std::string EvtVubHybrid::getName(){

  return "VUBHYBRID";     

}

EvtDecayBase* EvtVubHybrid::clone(){

  return new EvtVubHybrid;

}

void EvtVubHybrid::init(){

  // check that there are at least 3 arguments
  if (getNArg() < EvtVubHybrid::nParameters) {
    report(Severity::Error,"EvtVubHybrid") << "EvtVub generator expected "
				 << "at least " << EvtVubHybrid::nParameters
				 << " arguments but found: " << getNArg()
				 << "\nWill terminate execution!"<<endl;
    ::abort(); 
  } else if (getNArg() == EvtVubHybrid::nParameters) {
    report(Severity::Warning,"EvtVubHybrid") << "EvtVub: generate B -> Xu l nu events " 
				   << "without using the hybrid reweighting." 
				   << endl;
    _noHybrid = true;
  } else if (getNArg() < EvtVubHybrid::nParameters+EvtVubHybrid::nVariables) {
     report(Severity::Error,"EvtVubHybrid") << "EvtVub could not read number of bins for "
				  << "all variables used in the reweighting\n"
				  << "Will terminate execution!" << endl;
     ::abort();    
  }

  // check that there are 3 daughters
  checkNDaug(3);

  // read minimum required parameters from decay.dec
  _mb     = getArg(0);
  _a      = getArg(1);
  _alphas = getArg(2);

  // the maximum dGamma*p2 value depends on alpha_s only:
  const double dGMax0 = 3.;
  _dGMax = 0.21344+8.905*_alphas;
  if ( _dGMax < dGMax0 ) _dGMax = dGMax0;

  // for the Fermi Motion we need a B-Meson mass - but it's not critical
  // to get an exact value; in order to stay in the phase space for
  // B+- and B0 use the smaller mass
  
  static double mB0 = EvtPDL::getMaxMass(EvtPDL::getId("B0"));
  static double mBP = EvtPDL::getMaxMass(EvtPDL::getId("B+"));
  static double mB = (mB0<mBP?mB0:mBP);
  
  const double xlow = -_mb;
  const double xhigh = mB-_mb;
  const int aSize = 10000;

  EvtPFermi pFermi(_a,mB,_mb);
  // pf is the cumulative distribution normalized to 1.
  _pf.resize(aSize);
  for(int i=0;i<aSize;i++){
    double kplus = xlow + (double)(i+0.5)/((double)aSize)*(xhigh-xlow);
    if ( i== 0 )
      _pf[i] = pFermi.getFPFermi(kplus);
    else
      _pf[i] = _pf[i-1] + pFermi.getFPFermi(kplus);
  }
  for (size_t index=0; index<_pf.size(); index++) {
    _pf[index]/=_pf[_pf.size()-1];
  }

  _dGamma = new EvtVubdGamma(_alphas);
  
  if (_noHybrid) return;	// Without hybrid weighting, nothing else to do

  _nbins_mX = abs((int)getArg(3));
  _nbins_q2 = abs((int)getArg(4));
  _nbins_El = abs((int)getArg(5));

  int nextArg = EvtVubHybrid::nParameters + EvtVubHybrid::nVariables;

  _nbins = _nbins_mX*_nbins_q2*_nbins_El;	// Binning of weight table

  int expectArgs = nextArg + _nbins_mX +_nbins_q2 + _nbins_El + _nbins;

  if (getNArg() < expectArgs) {
    report(Severity::Error,"EvtVubHybrid")
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

void EvtVubHybrid::initProbMax(){
  noProbMax();
}

void EvtVubHybrid::decay( EvtParticle *p ){

  int j;
  // B+ -> u-bar specflav l+ nu
  
  EvtParticle *xuhad, *lepton, *neutrino;
  EvtVector4R p4;
  // R. Faccini 21/02/03
  // move the reweighting up , before also shooting the fermi distribution
  double x,z,p2;
  double sh=0.0;
  double mB,ml,xlow,xhigh,qplus;
  double El=0.0;
  double Eh=0.0;
  double kplus;
  double q2, mX;

  const double lp2epsilon=-10;
  bool rew(true);
  while(rew){
    
    p->initializePhaseSpace(getNDaug(),getDaugs());
    
    xuhad=p->getDaug(0);
    lepton=p->getDaug(1);
    neutrino=p->getDaug(2);
    
    mB = p->mass();
    ml = lepton->mass();
    
    xlow = -_mb;
    xhigh = mB-_mb;    
    
    // Fermi motion does not need to be computed inside the
    // tryit loop as m_b in Gamma0 does not need to be replaced by (m_b+kplus).
    // The difference however should be of the Order (lambda/m_b)^2 which is
    // beyond the considered orders in the paper anyway ...
    
    // for alpha_S = 0 and a mass cut on X_u not all values of kplus are 
    // possible. The maximum value is mB/2-_mb + sqrt(mB^2/4-_masscut^2)
    kplus = 2*xhigh;
    
    while( kplus >= xhigh || kplus <= xlow 
	   || (_alphas == 0 && kplus >= mB/2-_mb 
	       + sqrt(mB*mB/4-_masscut*_masscut))) {
      kplus = findPFermi(); //_pFermi->shoot();
      kplus = xlow + kplus*(xhigh-xlow);
    }
    qplus = mB-_mb-kplus;
    if( (mB-qplus)/2.<=ml) continue;
   
    int tryit = 1;
    while (tryit) {
      
      x = EvtRandom::Flat();
      z = EvtRandom::Flat(0,2);
      p2=EvtRandom::Flat();
      p2 = pow(10,lp2epsilon*p2);
      
      El = x*(mB-qplus)/2;
      if ( El > ml && El < mB/2) {
	
	Eh = z*(mB-qplus)/2+qplus;
	if ( Eh > 0 && Eh < mB ) {
	  
	  sh = p2*pow(mB-qplus,2)+2*qplus*(Eh-qplus)+qplus*qplus;
	  if ( sh > _masscut*_masscut
	       && mB*mB + sh - 2*mB*Eh > ml*ml) {
	    
	    double xran = EvtRandom::Flat();
	    
	    double y = _dGamma->getdGdxdzdp(x,z,p2)/_dGMax*p2;
	    
	    if ( y > 1 ) report(Severity::Warning,"EvtVubHybrid") <<"EvtVubHybrid decay probability > 1 found: " << y << endl;
 	    if ( y >= xran ) tryit = 0;
	  }
	}
      }
    }

    // compute all kinematic variables needed for reweighting (J. Dingfelder)
    mX = sqrt(sh);
    q2 = mB*mB + sh - 2*mB*Eh;

    // Reweighting in bins of mX, q2, El (J. Dingfelder)
    if (_nbins>0) {
      double xran1 = EvtRandom::Flat();
      double w = 1.0;
      if (!_noHybrid) w = getWeight(mX, q2, El); 
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
  ptmp = sqrt(Eh*Eh-sh);
  double pHB[4] = {Eh,ptmp*sttmp*cos(phH),ptmp*sttmp*sin(phH),ptmp*ctH};
  p4.set(pHB[0],pHB[1],pHB[2],pHB[3]);
  xuhad->init( getDaug(0), p4);

  if (_storeQplus ) {
    // cludge to store the hidden parameter q+ with the decay; 
    // the lifetime of the Xu is abused for this purpose.
    // tau = 1 ps corresponds to ctau = 0.3 mm -> in order to
    // stay well below BaBars sensitivity we take q+/(10000 GeV) which 
    // goes up to 0.0005 in the most extreme cases as ctau in mm.
    // To extract q+ back from the StdHepTrk its necessary to get
    // delta_ctau = Xu->anyDaughter->getVertexTime()-Xu->getVertexTime()
    // where these pseudo calls refere to the StdHep time stored at
    // the production vertex in the lab for each particle. The boost 
    // has to be reversed and the result is:
    //
    // q+ = delta_ctau * 10000 GeV/mm * Mass_Xu/Energy_Xu     
    //
    xuhad->setLifetime(qplus/10000.);
  }

  // calculate the W 4 vector in the B Meson restrframe

  double apWB = ptmp;
  double pWB[4] = {mB-Eh,-pHB[1],-pHB[2],-pHB[3]};

  // first go in the W restframe and calculate the lepton and
  // the neutrino in the W frame

  double mW2   = mB*mB + sh - 2*mB*Eh;
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
  // calculate the neutrino 4 vector in the W restframe
  //double pNW[4] = {sqrt(mW2)-pLW[0],-pLW[1],-pLW[2],-pLW[3]};
    
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


double EvtVubHybrid::findPFermi() {

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


double EvtVubHybrid::getWeight(double mX, double q2, double El) {

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


void EvtVubHybrid::readWeights(int startArg) {
  _weights  = new double[_nbins];

  double maxw = 0.0;
  for (int i = 0; i < _nbins; i++, startArg++) {
    _weights[i] = getArg(startArg);
    if (_weights[i] > maxw) maxw = _weights[i];
  }

  if (maxw == 0) {
    report(Severity::Error,"EvtVubHybrid") << "EvtVub generator expected at least one " 
				 << " weight > 0, but found none! " 
				 << "Will terminate execution!"<<endl;
    ::abort();
  }

  // rescale weights (to be in range 0..1)
  for (int i = 0; i < _nbins; i++) {
    _weights[i] /= maxw;
  }
}
