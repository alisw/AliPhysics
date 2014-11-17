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
// Module: EvtTSS.cc
//
// Description: Routine to decay tensor-> scalar scalar
//
// Modification history:
//
//    RYD     November 24, 1996         Module created
//
//------------------------------------------------------------------------
//
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtTSS.hh"
#include <string>

EvtTSS::~EvtTSS() {}

std::string EvtTSS::getName(){

  return "TSS";     

}


EvtDecayBase* EvtTSS::clone(){

  return new EvtTSS;

}

void EvtTSS::init(){

  // check that there are 0 arguments
  checkNArg(0);

  checkNDaug(2);

  checkSpinParent(EvtSpinType::TENSOR);

  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);


}

void EvtTSS::initProbMax() {

   setProbMax(1.0);

}      

void EvtTSS::decay( EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtVector4R moms1 = p->getDaug(0)->getP4(); 

  double norm = 1.0/(moms1.d3mag()*moms1.d3mag());

  vertex(0,norm*(p->epsTensor(0).cont1(EvtVector4C(moms1))*(moms1)));
  vertex(1,norm*(p->epsTensor(1).cont1(EvtVector4C(moms1))*(moms1)));
  vertex(2,norm*(p->epsTensor(2).cont1(EvtVector4C(moms1))*(moms1)));
  vertex(3,norm*(p->epsTensor(3).cont1(EvtVector4C(moms1))*(moms1)));
  vertex(4,norm*(p->epsTensor(4).cont1(EvtVector4C(moms1))*(moms1)));

  return ;
}
