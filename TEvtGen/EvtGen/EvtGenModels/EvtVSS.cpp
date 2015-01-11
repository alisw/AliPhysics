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
// Module: EvtVSS.cc
//
// Description: Routine to decay vector-> scalar scalar
//
// Modification history:
//
//    RYD       November 24, 1996       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtVSS.hh"
#include <string>

EvtVSS::~EvtVSS() {}

std::string EvtVSS::getName(){

  return "VSS"; 
    
}


EvtDecayBase* EvtVSS::clone(){

  return new EvtVSS;

}

void EvtVSS::init(){

  // check that there are 0 arguments
  checkNArg(0);

  // check that there are 2 daughters
  checkNDaug(2);

  // check the parent and daughter spins
  checkSpinParent(EvtSpinType::VECTOR);
  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);
}

void EvtVSS::initProbMax() {

   setProbMax(1.0);

}      

void EvtVSS::decay( EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtVector4R pDaug = p->getDaug(0)->getP4();
  
  double norm=1.0/pDaug.d3mag();

  for (int i=0; i<3; i++) vertex(i,norm*pDaug*(p->eps(i)));

  return;
}

