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
// Module: EvtVll.cc
//
// Description: The decay of a vector meson to two leptons,
//              or generally, two spin 1/2 particles.
//              E.g., J/psi -> e+ e-
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
#include "EvtGenModels/EvtVll.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

EvtVll::~EvtVll() {}

std::string EvtVll::getName(){

  return "VLL";     

}


EvtDecayBase* EvtVll::clone(){

  return new EvtVll;

}

void EvtVll::init(){

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::VECTOR);
  
  checkSpinDaughter(0,EvtSpinType::DIRAC);
  checkSpinDaughter(1,EvtSpinType::DIRAC);

}

void EvtVll::initProbMax(){

  setProbMax(1.0);

}

void EvtVll::decay(EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtParticle *l1, *l2;
  l1 = p->getDaug(0);
  l2 = p->getDaug(1);

  EvtVector4C l11, l12, l21, l22;
  l11=EvtLeptonVCurrent(l1->spParent(0),l2->spParent(0));
  l12=EvtLeptonVCurrent(l1->spParent(0),l2->spParent(1));
  l21=EvtLeptonVCurrent(l1->spParent(1),l2->spParent(0));
  l22=EvtLeptonVCurrent(l1->spParent(1),l2->spParent(1));

  EvtVector4C eps0=p->eps(0);
  EvtVector4C eps1=p->eps(1);
  EvtVector4C eps2=p->eps(2);

  double M2=p->mass();
  M2*=M2;
  double m2=l1->mass();
  m2*=m2;

  double norm=1.0/sqrt(2*M2+4*m2-4*m2*m2/M2);

  vertex(0,0,0,norm*(eps0*l11));
  vertex(0,0,1,norm*(eps0*l12));
  vertex(0,1,0,norm*(eps0*l21));
  vertex(0,1,1,norm*(eps0*l22));
  
  vertex(1,0,0,norm*(eps1*l11));
  vertex(1,0,1,norm*(eps1*l12));
  vertex(1,1,0,norm*(eps1*l21));
  vertex(1,1,1,norm*(eps1*l22));
  
  vertex(2,0,0,norm*(eps2*l11));
  vertex(2,0,1,norm*(eps2*l12));
  vertex(2,1,0,norm*(eps2*l21));
  vertex(2,1,1,norm*(eps2*l22));
  
  return;

}













