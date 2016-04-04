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
// Module: EvtHQET2.cc
//
// Description: Routine to implement semileptonic B->D*lnu & B->Dlnu 
//              decays according to the model HQET
//
//   Lange Nov9/01 adding Dlnu and possible (w-1)^2 term
//
//
// Modification history:
//
//    Marco Bomben     March 10, 2003       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <assert.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtHQET2.hh"
#include "EvtGenModels/EvtHQET2FF.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include <string>
using std::endl;

EvtHQET2::EvtHQET2():
  hqetffmodel(0)
  ,calcamp(0)
{}

EvtHQET2::~EvtHQET2() {
  delete hqetffmodel;
  hqetffmodel=0;
  delete calcamp;
  calcamp=0;
}

std::string EvtHQET2::getName(){

  return "HQET2";     

}



EvtDecayBase* EvtHQET2::clone(){

  return new EvtHQET2;

}


void EvtHQET2::decay( EvtParticle *p ){

  p->initializePhaseSpace(getNDaug(),getDaugs());
  calcamp->CalcAmp(p,_amp2,hqetffmodel);

}

void EvtHQET2::initProbMax(){

  EvtId parnum,mesnum,lnum,nunum;

  parnum = getParentId();
  mesnum = getDaug(0);
  lnum = getDaug(1);
  nunum = getDaug(2);

  double mymaxprob = calcamp->CalcMaxProb(parnum,mesnum,
                           lnum,nunum,hqetffmodel);

  setProbMax(mymaxprob);

}


void EvtHQET2::init(){

  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be X lepton neutrino
  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);

  EvtSpinType::spintype d1type = EvtPDL::getSpinType(getDaug(0));
  if ( d1type==EvtSpinType::SCALAR) {
    if ( getNArg()==2 ) {hqetffmodel = new EvtHQET2FF(getArg(0),getArg(1)); 
    calcamp = new EvtSemiLeptonicScalarAmp;} 
    else {
    report(Severity::Error,"EvtGen") << "HQET2 model for scalar meson daughters needs 2 arguments. Sorry."<<endl;
    ::abort();
  }  
  }
  else if ( d1type==EvtSpinType::VECTOR) {
    if ( getNArg()==4 ){ hqetffmodel = new EvtHQET2FF(getArg(0),getArg(1),getArg(2),getArg(3));
    calcamp = new EvtSemiLeptonicVectorAmp; }
    else  {
    report(Severity::Error,"EvtGen") << "HQET2 model for vector meson daughtersneeds 4 arguments. Sorry."<<endl;
    ::abort();
    }
  }
  else{
    report(Severity::Error,"EvtGen") << "HQET2 model handles only scalar and vector meson daughters. Sorry."<<endl;
    ::abort();
  }

  
}

