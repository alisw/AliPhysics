//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2002      INFN-Pisa
//
// Module: EvtVSSBMixCPT.cc
//
// Description:
//    Routine to decay vector-> scalar scalar with coherent BB-like mixing
//                              including CPT effects
//    Based on VSSBMIX
//
// Modification history:
//
//    F. Sandrelli, Fernando M-V March 03, 2002 
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenModels/EvtVSSBMixCPT.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include "EvtGenBase/EvtRandom.hh"
using std::endl;

EvtVSSBMixCPT::~EvtVSSBMixCPT() {}

std::string EvtVSSBMixCPT::getName(){
  return "VSS_BMIX";
}


EvtDecayBase* EvtVSSBMixCPT::clone(){
  return new EvtVSSBMixCPT;
}

void EvtVSSBMixCPT::init(){

  if ( getNArg()>4) checkNArg(14,12,8);

  if (getNArg()<1) {
    report(Severity::Error,"EvtGen") << "EvtVSSBMix generator expected "
                           << " at least 1 argument (deltam) but found:"<<getNArg()<<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }
  // check that we are asked to produced exactly 2 daughters
  //4 are allowed if they are aliased..
  checkNDaug(2,4);

  if ( getNDaug()==4) {
    if ( getDaug(0)!=getDaug(2)||getDaug(1)!=getDaug(3)){
      report(Severity::Error,"EvtGen") << "EvtVSSBMixCPT generator allows "
			     << " 4 daughters only if 1=3 and 2=4"
			     << " (but 3 and 4 are aliased "<<endl;
      report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
      ::abort();
    }
  }
  // check that we are asked to decay a vector particle into a pair
  // of scalar particles

  checkSpinParent(EvtSpinType::VECTOR);

  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);

  // check that our daughter particles are charge conjugates of each other
  if(!(EvtPDL::chargeConj(getDaug(0)) == getDaug(1))) {
    report(Severity::Error,"EvtGen") << "EvtVSSBMixCPT generator expected daughters "
			   << "to be charge conjugate." << endl
			   << "  Found " << EvtPDL::name(getDaug(0)).c_str() << " and "
			   << EvtPDL::name(getDaug(1)).c_str() << endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }
  // check that both daughter particles have the same lifetime
  if(EvtPDL::getctau(getDaug(0)) != EvtPDL::getctau(getDaug(1))) {
    report(Severity::Error,"EvtGen") << "EvtVSSBMixCPT generator expected daughters "
			   << "to have the same lifetime." << endl
			   << "  Found ctau = "
			   << EvtPDL::getctau(getDaug(0)) << " mm and "
			   << EvtPDL::getctau(getDaug(1)) << " mm" << endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }
  // precompute quantities that will be used to generate events
  // and print out a summary of parameters for this decay

  // mixing frequency in hbar/mm
  _freq= getArg(0)/EvtConst::c; 

  // deltaG 
  double gamma= 1/EvtPDL::getctau(getDaug(0));  // gamma/c (1/mm)
  _dGamma=0.0;
  double dgog=0.0;
  if ( getNArg() > 1 ) {
    dgog=getArg(1);
    _dGamma=dgog*gamma;
  }
  // q/p
  _qoverp = EvtComplex(1.0,0.0);
  if ( getNArg() > 2){
    _qoverp = EvtComplex(getArg(2),0.0); 
  } 
  if ( getNArg() > 3) {
    _qoverp = getArg(2)*EvtComplex(cos(getArg(3)),sin(getArg(3)));
  }
  _poverq=1.0/_qoverp;

  // decay amplitudes
  _A_f=EvtComplex(1.0,0.0);
  _Abar_f=EvtComplex(0.0,0.0);
  _A_fbar=_Abar_f;  // CPT conservation
  _Abar_fbar=_A_f;  // CPT conservation
  if ( getNArg() > 4){
    _A_f=getArg(4)*EvtComplex(cos(getArg(5)),sin(getArg(5)));    // this allows for DCSD
    _Abar_f=getArg(6)*EvtComplex(cos(getArg(7)),sin(getArg(7))); // this allows for DCSD 
    if ( getNArg() > 8 ){
      // CPT violation in decay 
      _A_fbar=getArg(8)*EvtComplex(cos(getArg(9)),sin(getArg(9)));
      _Abar_fbar=getArg(10)*EvtComplex(cos(getArg(11)),sin(getArg(11)));
    } else {
      // CPT conservation in decay
      _A_fbar=_Abar_f;
      _Abar_fbar=_A_f;
    }
  }

  // CPT violation in mixing
  _z = EvtComplex(0.0,0.0);
  if ( getNArg() > 12 ){
    _z = EvtComplex(getArg(12),getArg(13)); 
  }


  // some printout
  double tau= 1e12*EvtPDL::getctau(getDaug(0))/EvtConst::c; // in ps
  double dm= 1e-12*getArg(0); // B0/anti-B0 mass difference in hbar/ps
  double x= dm*tau; 
  double y= dgog*0.5; //y=dgamma/(2*gamma) 
  double qop2 = abs(_qoverp*_qoverp);
  _chib0_b0bar=qop2*(x*x+y*y)/(qop2*(x*x+y*y)+2+x*x-y*y);  // does not include CPT in mixing
  _chib0bar_b0=(1/qop2)*(x*x+y*y)/((1/qop2)*(x*x+y*y)+2+x*x-y*y); // does not include CPT in mixing 

  if ( verbose() ) {
    report(Severity::Info,"EvtGen") << "VSS_BMIXCPT will generate mixing and CPT/CP effects in mixing:"
			  << endl << endl
			  << "    " << EvtPDL::name(getParentId()).c_str() << " --> "
			  << EvtPDL::name(getDaug(0)).c_str() << " + "
			  << EvtPDL::name(getDaug(1)).c_str() << endl << endl
			  << "using parameters:" << endl << endl
			  << "  delta(m)  = " << dm << " hbar/ps" << endl
			  << "  _freq     = " << _freq << " hbar/mm" << endl
			  << "  dgog      = "  << dgog <<endl
			  << "  dGamma    = "  << _dGamma <<" hbar/mm" <<endl
			  << "  q/p       = " << _qoverp << endl  
			  << "  z         = " << _z << endl  
			  << "  tau       = " << tau << " ps" << endl
			  << "  x         = " << x << endl
			  << " chi(B0->B0bar) = " << _chib0_b0bar << endl
			  << " chi(B0bar->B0) = " << _chib0bar_b0 << endl 
			  << " Af         = " << _A_f << endl
			  << " Abarf      = " << _Abar_f << endl
			  << " Afbar      = " << _A_fbar << endl
			  << " Abarfbar   = " << _Abar_fbar << endl
			  << endl;
  }
}

void EvtVSSBMixCPT::initProbMax(){
  // this value is ok for reasonable values of all the parameters
  setProbMax(4.0);
}

void EvtVSSBMixCPT::decay( EvtParticle *p ){

  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  // generate a final state according to phase space

  double rndm= EvtRandom::random();

  if ( getNDaug()==4) {
    EvtId tempDaug[2];

    if ( rndm < 0.5 ) { tempDaug[0]=getDaug(0);  tempDaug[1]=getDaug(3); }
    else{ tempDaug[0]=getDaug(2);  tempDaug[1]=getDaug(1); }

    p->initializePhaseSpace(2,tempDaug);
  }
  else{ //nominal case.
    p->initializePhaseSpace(2,getDaugs());
  }

  EvtParticle *s1,*s2;

  s1 = p->getDaug(0);
  s2 = p->getDaug(1);
  //delete any daughters - if there are daughters, they
  //are from the initialization and will be redone later
  if ( s1->getNDaug() > 0 ) { s1->deleteDaughters();}
  if ( s2->getNDaug() > 0 ) { s2->deleteDaughters();}
  
  EvtVector4R p1= s1->getP4();
  EvtVector4R p2= s2->getP4();

  // throw a random number to decide if this final state should be mixed
  rndm= EvtRandom::random();
  int mixed= (rndm < 0.5) ? 1 : 0;

  // if this decay is mixed, choose one of the 2 possible final states
  // with equal probability (re-using the same random number)
  if(mixed==1) {
    EvtId mixedId= (rndm < 0.25) ? getDaug(0) : getDaug(1);
    EvtId mixedId2= mixedId;
    if (getNDaug()==4&&rndm<0.25)  mixedId2=getDaug(2);
    if (getNDaug()==4&&rndm>0.25)  mixedId2=getDaug(3);
    s1->init(mixedId, p1);
    s2->init(mixedId2, p2);
  }


  // if this decay is unmixed, choose one of the 2 possible final states
  // with equal probability (re-using the same random number)
  if(mixed==0) {
    EvtId unmixedId = (rndm < 0.75) ? getDaug(0) : getDaug(1);
    EvtId unmixedId2= (rndm < 0.75) ? getDaug(1) : getDaug(0);
    if (getNDaug()==4&&rndm<0.75)  unmixedId2=getDaug(3);
    if (getNDaug()==4&&rndm>0.75)  unmixedId2=getDaug(2);
    s1->init(unmixedId, p1);
    s2->init(unmixedId2, p2);
  }

  // choose a decay time for each final state particle using the
  // lifetime (which must be the same for both particles) in pdt.table
  // and calculate the lifetime difference for this event
  s1->setLifetime();
  s2->setLifetime();
  double dct= s1->getLifetime() - s2->getLifetime(); // in mm

  // Convention: _dGamma=GammaLight-GammaHeavy
  EvtComplex exp1(-0.25*_dGamma*dct,0.5*_freq*dct);

  /*
  //Find the flavor of the B that decayed first.
  EvtId firstDec = (dct > 0 ) ? s2->getId() : s1->getId();
 
  //This tags the flavor of the other particle at that time. 
  EvtId stateAtDeltaTeq0 = ( firstDec==B0 ) ? B0B : B0; 
  */
  EvtId stateAtDeltaTeq0 = (s2->getId()==B0) ? B0B : B0;

  // calculate the oscillation amplitude, based on wether this event is mixed or not
  EvtComplex osc_amp;

  //define some useful functions: (see BAD #188 eq. 39 for ref.) 
  EvtComplex gp=0.5*(exp(-1.0*exp1)+exp(exp1)); 
  EvtComplex gm=0.5*(exp(-1.0*exp1)-exp(exp1));
  EvtComplex sqz=sqrt(abs(1-_z*_z))*exp(EvtComplex(0,arg(1-_z*_z)/2));
  
  EvtComplex BB=gp+_z*gm;                // <B0|B0(t)> 
  EvtComplex barBB=-sqz*_qoverp*gm;      // <B0bar|B0(t)>
  EvtComplex BbarB=-sqz*_poverq*gm;      // <B0|B0bar(t)>
  EvtComplex barBbarB=gp-_z*gm;          // <B0bar|B0bar(t)>

  //
  if ( !mixed&&stateAtDeltaTeq0==B0 ) {
    osc_amp= BB*_A_f+barBB*_Abar_f;
  }
  if ( !mixed&&stateAtDeltaTeq0==B0B ) {
    osc_amp= barBbarB*_Abar_fbar+BbarB*_A_fbar;
  }

  if ( mixed&&stateAtDeltaTeq0==B0 ) {
    osc_amp=barBB*_Abar_fbar+BB*_A_fbar;  
  }
  if ( mixed&&stateAtDeltaTeq0==B0B ) {
    osc_amp=BbarB*_A_f+barBbarB*_Abar_f;
  }

  // store the amplitudes for each parent spin basis state
  double norm=1.0/p1.d3mag();
  vertex(0,norm*osc_amp*p1*(p->eps(0)));
  vertex(1,norm*osc_amp*p1*(p->eps(1)));
  vertex(2,norm*osc_amp*p1*(p->eps(2)));

  return ;
}

std::string EvtVSSBMixCPT::getParamName(int i) {
  switch(i) {
  case 0:
    return "deltaM";
  case 1:
    return "deltaGammaOverGamma";
  case 2:
    return "qOverP";
  case 3:
    return "qOverPPhase";
  case 4:
    return "Af";
  case 5:
    return "AfPhase";
  case 6:
    return "Abarf";
  case 7:
    return "AbarfPhase";
  case 8:
    return "Afbar";
  case 9:
    return "AfbarPhase";
  case 10:
    return "Abarfbar";
  case 11:
    return "AbarfbarPhase";
  case 12:
    return "Z";
  case 13:
    return "ZPhase";
  default:
    return "";
  }
}

std::string EvtVSSBMixCPT::getParamDefault(int i) {
  switch(i) {
  case 3:
    return "0.0";
  case 4:
    return "1.0";
  case 5:
    return "0.0";
  case 6:
    return "1.0";
  case 7:
    return "0.0";
  default:
    return "";
  }
}
