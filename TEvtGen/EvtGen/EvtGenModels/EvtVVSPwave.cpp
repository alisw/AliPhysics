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
// Module: EvtVVSPwave.cc
//
// Description: Routine to decay vector-> vector scalar in Partial-wave
//              Routine to decay a vector into a vector and scalar.  Started
//              by ryd on Aug 20, 1996.
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
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtVVSPwave.hh"
#include <string>
using std::endl;

EvtVVSPwave::~EvtVVSPwave() {}

std::string EvtVVSPwave::getName(){

  return "VVS_PWAVE";     

}


EvtDecayBase* EvtVVSPwave::clone(){

  return new EvtVVSPwave;

}

void EvtVVSPwave::init(){

  // check that there are 6 arguments
  checkNArg(6);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::VECTOR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);
}

void EvtVVSPwave::initProbMax() {

  //probmax is 1.0 for all possible decays I think!

   setProbMax(1.0);

}      

void EvtVVSPwave::decay( EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtComplex as(getArg(0)*cos(getArg(1)),getArg(0)*sin(getArg(1)));
  EvtComplex ap(getArg(2)*cos(getArg(3)),getArg(2)*sin(getArg(3)));
  EvtComplex ad(getArg(4)*cos(getArg(5)),getArg(4)*sin(getArg(5)));

  if (ap!=EvtComplex(0.0,0.0)) {
    report(Severity::Error,"EvtGen") << "In EvtVectorToVectorScalar.cc"<<endl;
    report(Severity::Error,"EvtGen") << "P wave not yet implemented!!"<<endl;
    ::abort();
  }
    
  EvtParticle *v;
  v=p->getDaug(0);

  EvtTensor4C d,g;
  
  g.setdiag(1.0,-1.0,-1.0,-1.0);

  d=ad*((1.0/(v->getP4().d3mag()*v->getP4().d3mag()))*
        EvtGenFunctions::directProd(v->getP4(),v->getP4())+(1/3.0)*g)+
    as*g;

  EvtVector4C ep0,ep1,ep2;  
  
  ep0=d.cont1(p->eps(0));
  ep1=d.cont1(p->eps(1));
  ep2=d.cont1(p->eps(2));

  vertex(0,0,ep0.cont(v->eps(0).conj()));
  vertex(0,1,ep0.cont(v->eps(1).conj()));
  vertex(0,2,ep0.cont(v->eps(2).conj()));
  
  vertex(1,0,ep1.cont(v->eps(0).conj()));
  vertex(1,1,ep1.cont(v->eps(1).conj()));
  vertex(1,2,ep1.cont(v->eps(2).conj()));
  
  vertex(2,0,ep2.cont(v->eps(0).conj()));
  vertex(2,1,ep2.cont(v->eps(1).conj()));
  vertex(2,2,ep2.cont(v->eps(2).conj()));

  return ;

}
