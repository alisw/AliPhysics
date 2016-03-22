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
// Module: EvtTauHadnu.cc
//
// Description: The leptonic decay of the tau meson.
//              E.g., tau- -> e- nueb nut
//
// Modification history:
//
//    RYD       January 17, 1997       Module created
//
//------------------------------------------------------------------------
//
#include <stdlib.h>
#include <iostream>
#include <string>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenModels/EvtTauHadnu.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtIdSet.hh"

using namespace std;

EvtTauHadnu::~EvtTauHadnu() {}

std::string EvtTauHadnu::getName(){

  return "TAUHADNU";     

}


EvtDecayBase* EvtTauHadnu::clone(){

  return new EvtTauHadnu;

}

void EvtTauHadnu::init() {

  // check that there are 0 arguments

  checkSpinParent(EvtSpinType::DIRAC);

  //the last daughter should be a neutrino
  checkSpinDaughter(getNDaug()-1,EvtSpinType::NEUTRINO);

  int i;
  for ( i=0; i<(getNDaug()-1);i++) {
    checkSpinDaughter(i,EvtSpinType::SCALAR);
  }

  bool validndaug=false;

  if ( getNDaug()==4 ) {
    //pipinu
    validndaug=true;
    checkNArg(7);
    _beta = getArg(0);
    _mRho = getArg(1);
    _gammaRho = getArg(2);
    _mRhopr = getArg(3);
    _gammaRhopr = getArg(4);
    _mA1 = getArg(5);
    _gammaA1 = getArg(6);
  }
  if ( getNDaug()==3 ) {
    //pipinu
    validndaug=true;
    checkNArg(5);
    _beta = getArg(0);
    _mRho = getArg(1);
    _gammaRho = getArg(2);
    _mRhopr = getArg(3);
    _gammaRhopr = getArg(4);
  }
  if ( getNDaug()==2 ) {
    //pipinu
    validndaug=true;
    checkNArg(0);
  }

  if ( !validndaug ) {
    report(Severity::Error,"EvtGen") << "Have not yet implemented this final state in TAUHADNUKS model" << endl;
    report(Severity::Error,"EvtGen") << "Ndaug="<<getNDaug() << endl;
    int id;
    for ( id=0; id<(getNDaug()-1); id++ ) 
      report(Severity::Error,"EvtGen") << "Daug " << id << " "<<EvtPDL::name(getDaug(id)).c_str() << endl;

  }

}

void EvtTauHadnu::initProbMax(){

  if ( getNDaug()==2 )  setProbMax(90.0);
  if ( getNDaug()==3 )  setProbMax(2500.0);
  if ( getNDaug()==4 )  setProbMax(30000.0);

}

void EvtTauHadnu::decay(EvtParticle *p){

  static EvtId TAUM=EvtPDL::getId("tau-");

  EvtIdSet thePis("pi+","pi-","pi0");
  EvtIdSet theKs("K+","K-");

  p->initializePhaseSpace(getNDaug(),getDaugs());
  
  EvtParticle *nut;
  nut = p->getDaug(getNDaug()-1);
  p->getDaug(0)->getP4();

  //get the leptonic current 
  EvtVector4C tau1, tau2;
  
  if (p->getId()==TAUM) {
    tau1=EvtLeptonVACurrent(nut->spParentNeutrino(),p->sp(0));
    tau2=EvtLeptonVACurrent(nut->spParentNeutrino(),p->sp(1));
  }
  else{
    tau1=EvtLeptonVACurrent(p->sp(0),nut->spParentNeutrino());
    tau2=EvtLeptonVACurrent(p->sp(1),nut->spParentNeutrino());
  }

  EvtVector4C hadCurr;
  bool foundHadCurr=false;
  if ( getNDaug() == 2 ) {
    hadCurr = p->getDaug(0)->getP4();
    foundHadCurr=true;
  }
  if ( getNDaug() == 3 ) {

    //pi pi0 nu with rho and rhopr resonance
    if ( thePis.contains(getDaug(0)) &&
	 thePis.contains(getDaug(1)) ) {

      EvtVector4R q1 = p->getDaug(0)->getP4();
      EvtVector4R q2 = p->getDaug(1)->getP4();
      double m1 = q1.mass();
      double m2 = q2.mass();
       
      hadCurr = Fpi( (q1+q2).mass2(), m1, m2 )  * (q1-q2);
      
      foundHadCurr = true;
    }

  }
  if ( getNDaug() == 4 ) {
    if ( thePis.contains(getDaug(0)) &&
	 thePis.contains(getDaug(1)) &&
	 thePis.contains(getDaug(2)) ) {
      //figure out which is the different charged pi
      //want it to be q3

      int diffPi(0),samePi1(0),samePi2(0);
      if ( getDaug(0) == getDaug(1) ) {diffPi=2; samePi1=0; samePi2=1;}
      if ( getDaug(0) == getDaug(2) ) {diffPi=1; samePi1=0; samePi2=2;}
      if ( getDaug(1) == getDaug(2) ) {diffPi=0; samePi1=1; samePi2=2;}

      EvtVector4R q1=p->getDaug(samePi1)->getP4();
      EvtVector4R q2=p->getDaug(samePi2)->getP4();
      EvtVector4R q3=p->getDaug(diffPi)->getP4();
      
      double m1 = q1.mass();
      double m2 = q2.mass();
      double m3 = q3.mass();
      
      EvtVector4R Q = q1 + q2 + q3;
      double Q2 = Q.mass2();
      double _mA12 = _mA1*_mA1;

      double _gammaA1X = _gammaA1 * gFunc( Q2, samePi1 )/gFunc( _mA12, samePi1 );

      EvtComplex denBW_A1( _mA12 - Q2, -1.*_mA1*_gammaA1X );
      EvtComplex BW_A1 = _mA12 / denBW_A1;

      hadCurr = BW_A1 * ( ((q1-q3)-(Q*(Q*(q1-q3))/Q2)) * Fpi( (q1+q3).mass2(), m1, m3) + 
			  ((q2-q3)-(Q*(Q*(q2-q3))/Q2)) * Fpi( (q2+q3).mass2(), m2, m3) ); 

      foundHadCurr = true;
      
    }


  }



  if ( !foundHadCurr ) {
    report(Severity::Error,"EvtGen") << "Have not yet implemented this final state in TAUHADNUKS model" << endl;
    report(Severity::Error,"EvtGen") << "Ndaug="<<getNDaug() << endl;
    int id;
    for ( id=0; id<(getNDaug()-1); id++ ) 
      report(Severity::Error,"EvtGen") << "Daug " << id << " "<<EvtPDL::name(getDaug(id)).c_str() << endl;

  }

  
  vertex(0,tau1*hadCurr);
  vertex(1,tau2*hadCurr);
  

  
  return;

}

double EvtTauHadnu::gFunc(double Q2, int dupD) {
  
  double mpi= EvtPDL::getMeanMass(getDaug(dupD));
  double mpi2 = pow( mpi,2.);
  if ( Q2 < pow(_mRho + mpi, 2.) ) {
    double arg = Q2-9.*mpi2;
    return 4.1 * pow(arg,3.) * (1. - 3.3*arg + 5.8*pow(arg,2.));
  }
  else 
    return Q2 * (1.623 + 10.38/Q2 - 9.32/pow(Q2,2.) + 0.65/pow(Q2,3.));
}

EvtComplex EvtTauHadnu::Fpi( double s, double xm1, double xm2 ) {

  EvtComplex BW_rho = BW( s, _mRho, _gammaRho, xm1, xm2 );
  EvtComplex BW_rhopr = BW( s, _mRhopr, _gammaRhopr, xm1, xm2 );
  
  
  return (BW_rho + _beta*BW_rhopr)/(1.+_beta);
}

EvtComplex EvtTauHadnu::BW( double s, double m, double gamma, double xm1, double xm2 ) {
  
  double m2 = pow( m, 2.);
  
  if ( s > pow( xm1+xm2, 2.) ) {
    double qs = sqrt( fabs( (s-pow(xm1+xm2,2.)) * (s-pow(xm1-xm2,2.)) ) ) / sqrt(s); 
    double qm = sqrt( fabs( (m2-pow(xm1+xm2,2.)) * (m2-pow(xm1-xm2,2.)) ) ) / m;
    
    gamma *= m2/s * pow( qs/qm, 3.); 
  }
  else
    gamma = 0.;

  EvtComplex denBW( m2 - s, -1.* sqrt(s) * gamma );
  
  
  return m2 / denBW;
}





