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
// Module: EvtSLN.cc
//
// Description: B ==> tau + nu
//
// Modification history:
//
//    RYD/SHY   April 23, 1997           Module created
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
#include "EvtGenModels/EvtSLN.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

EvtSLN::~EvtSLN() {}

std::string EvtSLN::getName(){

  return "SLN";     

}

EvtDecayBase* EvtSLN::clone(){

  return new EvtSLN;

}


void EvtSLN::init(){

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(2);
    
  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::DIRAC);
  checkSpinDaughter(1,EvtSpinType::NEUTRINO);

}


void EvtSLN::initProbMax(){

  double M=EvtPDL::getMeanMass(getParentId());
  double m=EvtPDL::getMeanMass(getDaug(0));

  double probMax=8.0*(M*M-m*m)*m*m;

  setProbMax(probMax);

}


void EvtSLN::decay(EvtParticle *p){

  static EvtId EM=EvtPDL::getId("e-");
  static EvtId MUM=EvtPDL::getId("mu-");
  static EvtId TAUM=EvtPDL::getId("tau-");

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtParticle *l, *nul;
  l= p->getDaug(0);
  nul= p->getDaug(1);

  EvtVector4R p4_p;
  p4_p.set(p->mass(),0.0,0.0,0.0);

  EvtVector4C l1, l2;
  
  if (getDaug(0)==TAUM || getDaug(0)==MUM || getDaug(0)==EM) {
    l1=EvtLeptonVACurrent(l->spParent(0),nul->spParentNeutrino());
    l2=EvtLeptonVACurrent(l->spParent(1),nul->spParentNeutrino());
  }
  else{
    l1=EvtLeptonVACurrent(nul->spParentNeutrino(),l->spParent(0));
    l2=EvtLeptonVACurrent(nul->spParentNeutrino(),l->spParent(1));
  }

  vertex(0,p4_p*l1);
  vertex(1,p4_p*l2);
  
  return;

}


