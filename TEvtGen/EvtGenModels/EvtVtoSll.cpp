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
// Description: The decay of a vector meson to a scalar and a 
//              lepton pair. E.g. D_s*+ -> D_s+ e+ e-
//           
//
// Modification history:
//
//    RYD       February 28, 2009       Module created
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
#include "EvtGenModels/EvtVtoSll.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"

EvtVtoSll::~EvtVtoSll() {}

std::string EvtVtoSll::getName(){

  return "VTOSLL";     

}


EvtDecayBase* EvtVtoSll::clone(){

  return new EvtVtoSll;

}

void EvtVtoSll::init(){

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(3);

  checkSpinParent(EvtSpinType::VECTOR);
  
  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::DIRAC);

}

void EvtVtoSll::initProbMax(){

  //setProbMax(1.0);

}

void EvtVtoSll::decay(EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtParticle *l1, *l2;
  l1 = p->getDaug(1);
  l2 = p->getDaug(2);

  EvtVector4C l11, l12, l21, l22;
  l11=EvtLeptonVCurrent(l1->spParent(0),l2->spParent(0));
  l12=EvtLeptonVCurrent(l1->spParent(0),l2->spParent(1));
  l21=EvtLeptonVCurrent(l1->spParent(1),l2->spParent(0));
  l22=EvtLeptonVCurrent(l1->spParent(1),l2->spParent(1));

  EvtVector4C eps0=p->eps(0);
  EvtVector4C eps1=p->eps(1);
  EvtVector4C eps2=p->eps(2);

  EvtVector4R P=p->getP4Restframe();
  EvtVector4R k=l1->getP4()+l2->getP4();
  double k2=k*k;

  EvtTensor4C T(dual(EvtGenFunctions::directProd(P,(1.0/k2)*k)));

  double M2=p->mass();
  M2*=M2;
  double m2=l1->mass();
  m2*=m2;

  double norm=1.0/sqrt(2*M2+4*m2-4*m2*m2/M2);

  vertex(0,0,0,norm*(eps0*T.cont2(l11)));
  vertex(0,0,1,norm*(eps0*T.cont2(l12)));
  vertex(0,1,0,norm*(eps0*T.cont2(l21)));
  vertex(0,1,1,norm*(eps0*T.cont2(l22)));
  
  vertex(1,0,0,norm*(eps1*T.cont2(l11)));
  vertex(1,0,1,norm*(eps1*T.cont2(l12)));
  vertex(1,1,0,norm*(eps1*T.cont2(l21)));
  vertex(1,1,1,norm*(eps1*T.cont2(l22)));
  
  vertex(2,0,0,norm*(eps2*T.cont2(l11)));
  vertex(2,0,1,norm*(eps2*T.cont2(l12)));
  vertex(2,1,0,norm*(eps2*T.cont2(l21)));
  vertex(2,1,1,norm*(eps2*T.cont2(l22)));
  
  return;

}













