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
// Module: EvtSLPole.cc
//
// Description: Routine to implement semileptonic decays according
//              to light cone sum rules
//
// Modification history:
//
//    DJL       April 23, 1998       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtSLPole.hh"
#include "EvtGenModels/EvtSLPoleFF.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicTensorAmp.hh"
#include <string>

EvtSLPole::~EvtSLPole() {}

std::string EvtSLPole::getName(){

  return "SLPOLE";     

}


EvtDecayBase* EvtSLPole::clone(){

  return new EvtSLPole;

}

void EvtSLPole::decay( EvtParticle *p ){

  p->initializePhaseSpace(getNDaug(),getDaugs(),_resetDaughterTree);
  calcamp->CalcAmp(p,_amp2,SLPoleffmodel);
  return;
}

void EvtSLPole::initProbMax(){

EvtId parnum,mesnum,lnum,nunum;

parnum = getParentId();
mesnum = getDaug(0);
lnum = getDaug(1);
nunum = getDaug(2);

double mymaxprob = calcamp->CalcMaxProb(parnum,mesnum,
                           lnum,nunum,SLPoleffmodel);

setProbMax(mymaxprob);

}


void EvtSLPole::init(){
  
  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be X lepton neutrino

  checkSpinParent(EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);

  EvtSpinType::spintype mesontype=EvtPDL::getSpinType(getDaug(0));

  SLPoleffmodel = new EvtSLPoleFF(getNArg(),getArgs());
  
  if ( mesontype==EvtSpinType::SCALAR ) { 
    calcamp = new EvtSemiLeptonicScalarAmp; 
  }
  if ( mesontype==EvtSpinType::VECTOR ) { 
    calcamp = new EvtSemiLeptonicVectorAmp; 
  }
  if ( mesontype==EvtSpinType::TENSOR ) { 
    calcamp = new EvtSemiLeptonicTensorAmp; 
  }

  _resetDaughterTree=false;
  if ( getArgStr(getNArg()-1) == "true") _resetDaughterTree=true;
  
}

