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
// Module: EvtSVVHelAmp.cc
//
// Description: Routine to decay scalar -> 2 vectors
//              by specifying the helicity amplitudes
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
#include "EvtGenBase/EvtVector3C.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "EvtGenBase/EvtTensor3C.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtSVVHelAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>

EvtSVVHelAmp::~EvtSVVHelAmp() {}

std::string EvtSVVHelAmp::getName(){

  return "SVV_HELAMP";     

}


EvtDecayBase* EvtSVVHelAmp::clone(){

  return new EvtSVVHelAmp;

}

void EvtSVVHelAmp::init(){

  // check that there are 6 arguments
  checkNArg(6);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::VECTOR);

}


void EvtSVVHelAmp::initProbMax(){

  setProbMax(getArg(0)*getArg(0)+getArg(2)*getArg(2)+getArg(4)*getArg(4));

}


void EvtSVVHelAmp::decay( EvtParticle *p){

  SVVHel(p,_amp2,getDaug(0),getDaug(1),
	 EvtComplex(getArg(0)*cos(getArg(1)),getArg(0)*sin(getArg(1))),
	 EvtComplex(getArg(2)*cos(getArg(3)),getArg(2)*sin(getArg(3))),
	 EvtComplex(getArg(4)*cos(getArg(5)),getArg(4)*sin(getArg(5))));
				 
  return ;

}


void EvtSVVHelAmp::SVVHel(EvtParticle *parent,EvtAmp& amp,EvtId n_v1,EvtId n_v2,
               const EvtComplex& hp,const EvtComplex& h0,
               const EvtComplex& hm){

  //  Routine to decay a vector into a vector and scalar.  Started
  //  by ryd on Oct 17, 1996.
    
  int tndaug = 2;
  EvtId tdaug[2];
  tdaug[0] = n_v1;
  tdaug[1] = n_v2;


  parent->initializePhaseSpace(tndaug,tdaug);

  EvtParticle *v1,*v2;
  v1 = parent->getDaug(0);
  v2 = parent->getDaug(1);

  EvtVector4R momv1 = v1->getP4();
  //EvtVector4R momv2 = v2->getP4();

  EvtVector3R v1dir(momv1.get(1),momv1.get(2),momv1.get(3));
  v1dir=v1dir/v1dir.d3mag();

  EvtComplex a=-0.5*(hp+hm);
  EvtComplex b=EvtComplex(0.0,0.5)*(hp-hm);
  EvtComplex c=h0+0.5*(hp+hm);

  EvtTensor3C M=a*EvtTensor3C::id()+
    b*EvtGenFunctions::eps(v1dir)+
    c*EvtGenFunctions::directProd(v1dir,v1dir);

  EvtVector3C t0=M.cont1(v1->eps(0).vec().conj());
  EvtVector3C t1=M.cont1(v1->eps(1).vec().conj());
  EvtVector3C t2=M.cont1(v1->eps(2).vec().conj());

  EvtVector3C eps0=v2->eps(0).vec().conj();
  EvtVector3C eps1=v2->eps(1).vec().conj();
  EvtVector3C eps2=v2->eps(2).vec().conj();

  amp.vertex(0,0,t0*eps0);
  amp.vertex(0,1,t0*eps1);
  amp.vertex(0,2,t0*eps2);

  amp.vertex(1,0,t1*eps0);
  amp.vertex(1,1,t1*eps1);
  amp.vertex(1,2,t1*eps2);

  amp.vertex(2,0,t2*eps0);
  amp.vertex(2,1,t2*eps1);
  amp.vertex(2,2,t2*eps2);

  return ;

}

std::string EvtSVVHelAmp::getParamName(int i) {
  switch(i) {
  case 0:
    return "plusHelAmp";
  case 1:
    return "plusHelAmpPhase";
  case 2:
    return "zeroHelAmp";
  case 3:
    return "zeroHelAmpPhase";
  case 4:
    return "minusHelAmp";
  case 5:
    return "minusHelAmpPhase";
  default:
    return "";
  }
}

std::string EvtSVVHelAmp::getParamDefault(int i) {
  switch(i) {
  case 0:
    return "1.0";
  case 1:
    return "0.0";
  case 2:
    return "1.0";
  case 3:
    return "0.0";
  case 4:
    return "1.0";
  case 5:
    return "0.0";
  default:
    return "";
  }
}
