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
// Module: EvtGenModels/EvtBToVlnuBall.cc
//
// Description:   B->Xu l nu with the Ball/Zwicky decay model
//                Xu is a vector (rho, rho0, omega)
//
//
// Modification history:
//
//    Wells Wulsin      2008 Aug 14         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <assert.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBToVlnuBall.hh"
#include "EvtGenModels/EvtBToVlnuBallFF.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include <string>
using std::endl;

EvtBToVlnuBall::EvtBToVlnuBall():
  _Ballmodel(0)
  ,_calcamp(0)
{}

EvtBToVlnuBall::~EvtBToVlnuBall() {
  delete _Ballmodel;
  _Ballmodel=0;
  delete _calcamp;
  _calcamp=0;
}


std::string EvtBToVlnuBall::getName(){
  return "BTOVLNUBALL";
}


EvtDecayBase* EvtBToVlnuBall::clone(){
  return new EvtBToVlnuBall;
}


void EvtBToVlnuBall::decay( EvtParticle *p ){

  p->initializePhaseSpace(getNDaug(),getDaugs());
  _calcamp->CalcAmp(p,_amp2,_Ballmodel);




}

void EvtBToVlnuBall::initProbMax(){

EvtId parnum,mesnum,lnum,nunum;

parnum = getParentId();
mesnum = getDaug(0);
lnum = getDaug(1);
nunum = getDaug(2);

double mymaxprob = _calcamp->CalcMaxProb(parnum,mesnum,
                           lnum,nunum,_Ballmodel);

setProbMax(mymaxprob);

}


void EvtBToVlnuBall::init(){

  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be X lepton neutrino
  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);

  EvtSpinType::spintype d1type = EvtPDL::getSpinType(getDaug(0));
  if ( d1type==EvtSpinType::VECTOR) {
    checkNArg(8);  // the number of arguments needed for the Ball model
    _Ballmodel = new EvtBToVlnuBallFF(getArg(0),getArg(1),getArg(2),getArg(3),getArg(4),getArg(5),getArg(6),getArg(7));
    _calcamp = new EvtSemiLeptonicVectorAmp; 
  } else {
    report(Severity::Error,"EvtGen") << "Ball model handles only vector meson daughters. Sorry." << endl;
    ::abort();
  }

  
}

