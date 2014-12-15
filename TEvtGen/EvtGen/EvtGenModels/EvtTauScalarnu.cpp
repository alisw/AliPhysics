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
// Module: EvtTauScalarnu.cc
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
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <iostream>
#include <string>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenModels/EvtTauScalarnu.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

EvtTauScalarnu::~EvtTauScalarnu() {}

std::string EvtTauScalarnu::getName(){

  return "TAUSCALARNU";     

}


EvtDecayBase* EvtTauScalarnu::clone(){

  return new EvtTauScalarnu;

}

void EvtTauScalarnu::init(){

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::DIRAC);

  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::NEUTRINO);

}

void EvtTauScalarnu::initProbMax(){

  setProbMax(90.0);

}

void EvtTauScalarnu::decay(EvtParticle *p){

  static EvtId TAUM=EvtPDL::getId("tau-");
  p->initializePhaseSpace(getNDaug(),getDaugs());
  
  EvtParticle *nut;
  nut = p->getDaug(1);
  EvtVector4R momscalar = p->getDaug(0)->getP4();
 
  EvtVector4C tau1, tau2;
  
  if (p->getId()==TAUM) {
    tau1=EvtLeptonVACurrent(nut->spParentNeutrino(),p->sp(0));
    tau2=EvtLeptonVACurrent(nut->spParentNeutrino(),p->sp(1));
  }
  else{
    tau1=EvtLeptonVACurrent(p->sp(0),nut->spParentNeutrino());
    tau2=EvtLeptonVACurrent(p->sp(1),nut->spParentNeutrino());
  }
  
  vertex(0,tau1*momscalar);
  vertex(1,tau2*momscalar);

  return;

}

