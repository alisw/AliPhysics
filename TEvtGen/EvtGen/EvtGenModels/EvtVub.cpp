//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtVub.cc
//
// Description: Routine to decay a particle according th phase space 
//
// Modification history:
//
//    Sven Menke       January 17, 2001       Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtVub.hh"
#include <string>
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenModels/EvtPFermi.hh"
#include "EvtGenModels/EvtVubdGamma.hh"
#include "EvtGenBase/EvtRandom.hh"
using std::endl;

EvtVub::~EvtVub() {
  if (_dGamma) delete _dGamma;
  if (_masses) delete [] _masses;
  if (_weights) delete [] _weights;
}

std::string EvtVub::getName(){

  return "VUB";     

}

EvtDecayBase* EvtVub::clone(){

  return new EvtVub;

}


void EvtVub::init(){

  // check that there are at least 6 arguments

  if (getNArg()<6) {

    report(Severity::Error,"EvtGen") << "EvtVub generator expected "
                           << " at least 6 arguments (mb,a,alpha_s,Nbins,m1,w1,...) but found: "
			   <<getNArg()<<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();

  }
  
  _mb      	= getArg(0);
  _a       	= getArg(1);
  _alphas  	= getArg(2);
  _nbins        = abs((int)getArg(3));
  _storeQplus   = (getArg(3)<0?1:0);
  _masses       = new double[_nbins];
  _weights      = new double[_nbins];
 
  if (getNArg()-4 != 2*_nbins) {
    report(Severity::Error,"EvtGen") << "EvtVub generator expected " 
                           << _nbins << " masses and weights but found: "
			   <<(getNArg()-4)/2 <<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }
  int i,j = 4;
  double maxw = 0.;
  for (i=0;i<_nbins;i++) {
    _masses[i] = getArg(j++);
    if (i>0 && _masses[i] <= _masses[i-1]) {
      report(Severity::Error,"EvtGen") << "EvtVub generator expected " 
			     << " mass bins in ascending order!"
			     << "Will terminate execution!"<<endl;
      ::abort();
    }
    _weights[i] = getArg(j++);
    if (_weights[i] < 0) {
      report(Severity::Error,"EvtGen") << "EvtVub generator expected " 
			     << " weights >= 0, but found: " 
			     <<_weights[i] <<endl;
      report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
      ::abort();
    }
    if ( _weights[i] > maxw ) maxw = _weights[i];
  }
  if (maxw == 0) {
    report(Severity::Error,"EvtGen") << "EvtVub generator expected at least one " 
			   << " weight > 0, but found none! " 
			   << "Will terminate execution!"<<endl;
    ::abort();
  }
  for (i=0;i<_nbins;i++) _weights[i]/=maxw;

  // the maximum dGamma*p2 value depends on alpha_s only:

  const double dGMax0 = 3.;
  _dGMax = 0.21344+8.905*_alphas;
  if ( _dGMax < dGMax0 ) _dGMax = dGMax0;

  // for the Fermi Motion we need a B-Meson mass - but it's not critical
  // to get an exact value; in order to stay in the phase space for
  // B+- and B0 use the smaller mass

  EvtId BP=EvtPDL::getId("B+");
  EvtId B0=EvtPDL::getId("B0");
  
  double mB0 = EvtPDL::getMaxMass(B0);
  double mBP = EvtPDL::getMaxMass(BP);

  double mB = (mB0<mBP?mB0:mBP);
  
  const double xlow = -_mb;
  const double xhigh = mB-_mb;
  const int aSize = 10000;

  EvtPFermi pFermi(_a,mB,_mb);
  // pf is the cumulative distribution
  // normalized to 1.
  _pf.resize(aSize);
  for(i=0;i<aSize;i++){
    double kplus = xlow + (double)(i+0.5)/((double)aSize)*(xhigh-xlow);
    if ( i== 0 )
      _pf[i] = pFermi.getFPFermi(kplus);
    else
      _pf[i] = _pf[i-1] + pFermi.getFPFermi(kplus);
  }
  for (size_t index=0; index<_pf.size(); index++) {
    _pf[index]/=_pf[_pf.size()-1];
  }

  //  static EvtHepRandomEngine myEngine;

  //  _pFermi = new RandGeneral(myEngine,pf,aSize,0);
  _dGamma = new EvtVubdGamma(_alphas);
  
  // check that there are 3 daughters
  checkNDaug(3);
}

void EvtVub::initProbMax(){

  noProbMax();

}

void EvtVub::decay( EvtParticle *p ){

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
    // possible. The maximum value is mB/2-_mb + sqrt(mB^2/4-_masses[0]^2)
    kplus = 2*xhigh;
    
    while( kplus >= xhigh || kplus <= xlow 
	   || (_alphas == 0 && kplus >= mB/2-_mb 
	       + sqrt(mB*mB/4-_masses[0]*_masses[0]))) {
      kplus = findPFermi(); //_pFermi->shoot();
      kplus = xlow + kplus*(xhigh-xlow);
    }
    qplus = mB-_mb-kplus;
    if( (mB-qplus)/2.<=ml)continue;
   
    int tryit = 1;
    while (tryit) {
      
      x = EvtRandom::Flat();
      z = EvtRandom::Flat(0,2);
      p2=EvtRandom::Flat();
      p2 = pow(10.0,lp2epsilon*p2);
      
      El = x*(mB-qplus)/2;
      if ( El > ml && El < mB/2) {
	
	Eh = z*(mB-qplus)/2+qplus;
	if ( Eh > 0 && Eh < mB ) {
	  
	  sh = p2*pow(mB-qplus,2)+2*qplus*(Eh-qplus)+qplus*qplus;
	  if ( sh > _masses[0]*_masses[0]
	       && mB*mB + sh - 2*mB*Eh > ml*ml) {
	    
	    double xran = EvtRandom::Flat();
	    
	    double y = _dGamma->getdGdxdzdp(x,z,p2)/_dGMax*p2;
	    
	    if ( y > 1 ) report(Severity::Warning,"EvtGen")<<"EvtVub decay probability > 1 found: " << y << endl;
	    if ( y >= xran ) tryit = 0;
	  }
	}
      }
    }
    // reweight the Mx distribution
    if(_nbins>0){
      double xran1 = EvtRandom::Flat();
      double m = sqrt(sh);j=0;
      while ( j < _nbins && m > _masses[j] ) j++; 
      double w = _weights[j-1]; 
      if ( w >= xran1 ) rew = false;
    } else {
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
  double phH = EvtRandom::Flat(0,2*EvtConst::pi);
  double phL = EvtRandom::Flat(0,2*EvtConst::pi);

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


double EvtVub::findPFermi() {

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
