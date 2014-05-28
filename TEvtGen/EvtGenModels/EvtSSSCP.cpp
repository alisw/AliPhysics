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
// Module: EvtSSSCP.cc
//
// Description: Routine to decay scalar -> 2 scalars
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
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtSSSCP.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"

EvtSSSCP::~EvtSSSCP() {}

std::string EvtSSSCP::getName(){

  return "SSS_CP";     

}


EvtDecayBase* EvtSSSCP::clone(){

  return new EvtSSSCP;

}

void EvtSSSCP::init(){

  // check that there are 7 arguments
  checkNArg(7);
  checkNDaug(2);
  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);

}

void EvtSSSCP::initProbMax(){

  //This is probably not quite right, but it should do as a start...
  //Anders

  setProbMax(2*(getArg(3)*getArg(3)+getArg(5)*getArg(5)));

}

void EvtSSSCP::decay( EvtParticle *p ){

  //added by Lange Jan4,2000
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  double t;
  EvtId other_b;

  EvtCPUtil::getInstance()->OtherB(p,t,other_b,0.5);

  p->initializePhaseSpace(getNDaug(),getDaugs());


  EvtComplex amp;

  EvtComplex A,Abar;
  
  A=EvtComplex(getArg(3)*cos(getArg(4)),getArg(3)*sin(getArg(4)));
  Abar=EvtComplex(getArg(5)*cos(getArg(6)),getArg(5)*sin(getArg(6)));
   
  if (other_b==B0B){
    amp=A*cos(getArg(1)*t/(2*EvtConst::c))+
      EvtComplex(cos(-2.0*getArg(0)),sin(-2.0*getArg(0)))*
      getArg(2)*EvtComplex(0.0,1.0)*Abar*sin(getArg(1)*t/(2*EvtConst::c));
  }
  if (other_b==B0){
    amp=A*EvtComplex(cos(2.0*getArg(0)),sin(2.0*getArg(0)))*
      EvtComplex(0.0,1.0)*sin(getArg(1)*t/(2*EvtConst::c))+       
      getArg(2)*Abar*cos(getArg(1)*t/(2*EvtConst::c));
  }
  
  vertex(amp);
  
  return ;
}

std::string EvtSSSCP::getParamName(int i) {
  switch(i) {
  case 0:
    return "weakPhase";
  case 1:
    return "deltaM";
  case 2:
    return "finalStateCP";
  case 3:
    return "Af";
  case 4:
    return "AfPhase";
  case 5:
    return "Abarf";
  case 6:
    return "AbarfPhase";
  default:
    return "";
  }
}

std::string EvtSSSCP::getParamDefault(int i) {
  switch(i) {
  case 3:
    return "1.0";
  case 4:
    return "0.0";
  case 5:
    return "1.0";
  case 6:
    return "0.0";
  default:
    return "";
  }
}
