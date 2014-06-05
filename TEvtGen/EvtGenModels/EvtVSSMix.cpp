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
// Module: EvtVSSMix.cc
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
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenModels/EvtVSSMix.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"

EvtVSSMix::~EvtVSSMix() {}

std::string EvtVSSMix::getName(){

  return "VSS_MIX";     

}


EvtDecayBase* EvtVSSMix::clone(){

  return new EvtVSSMix;

}

void EvtVSSMix::init(){

  // check that there are 1 arguments
  checkNArg(1);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::VECTOR);
    
  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);

}

void EvtVSSMix::initProbMax(){

  setProbMax(0.5);

}

void EvtVSSMix::decay( EvtParticle *p ){

  //added by Lange Jan4,2000
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  p->initializePhaseSpace(getNDaug(),getDaugs());
  EvtParticle *s1,*s2;
  s1 = p->getDaug(0);
  s2 = p->getDaug(1);
  EvtVector4R s1mom = s1->getP4();

  double t1,t2,dm;

  s1->setLifetime();
  s2->setLifetime();

  t1=s1->getLifetime();
  t2=s2->getLifetime();

  //dm should probably be a parameter to this model.

  dm=getArg(0)/EvtConst::c;

  EvtId d1,d2;

  d1=s1->getId();
  d2=s2->getId();

  double mix_amp=0.;
  if (d1==B0&&d2==B0B) mix_amp=cos(0.5*dm*(t1-t2));
  if (d1==B0B&&d2==B0) mix_amp=cos(0.5*dm*(t1-t2));
  if (d1==B0&&d2==B0) mix_amp=sin(0.5*dm*(t1-t2));
  if (d1==B0B&&d2==B0B) mix_amp=sin(0.5*dm*(t1-t2));

  double norm=1.0/s1mom.d3mag();

  vertex(0,norm*mix_amp*s1mom*(p->eps(0)));
  vertex(1,norm*mix_amp*s1mom*(p->eps(1)));
  vertex(2,norm*mix_amp*s1mom*(p->eps(2)));

  return ;
}

std::string EvtVSSMix::getParamName(int i) {
  switch(i) {
  case 0:
    return "deltaM";
  default:
    return "";
  }
}
