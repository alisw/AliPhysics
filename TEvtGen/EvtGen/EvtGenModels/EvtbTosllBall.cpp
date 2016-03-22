//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2000      Caltech, UCSB
//
// Module: EvtbTosllBall.cc
//
// Description: Routine to implement b->sll decays according to Ball et al. 
//
// Modification history:
//
//    Ryd     January 5, 2000        Module created
//
//    jjhollar October 7, 2005       Option to select form factors at runtime
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtbTosllBall.hh"
#include "EvtGenModels/EvtbTosllBallFF.hh"
#include "EvtGenModels/EvtbTosllAmp.hh"
#include "EvtGenModels/EvtbTosllScalarAmp.hh"
#include "EvtGenModels/EvtbTosllVectorAmp.hh"

#include <string>
using std::endl;

EvtbTosllBall::~EvtbTosllBall() {
  delete _calcamp;
  delete _ballffmodel;
}

std::string EvtbTosllBall::getName(){

  return "BTOSLLBALL";     
}


EvtDecayBase* EvtbTosllBall::clone(){

  return new EvtbTosllBall;

}

void EvtbTosllBall::decay( EvtParticle *p ){

  setWeight(p->initializePhaseSpace(getNDaug(),getDaugs(),false,
                                    _poleSize,1,2));

  _calcamp->CalcAmp(p,_amp2,_ballffmodel);
  
}


void EvtbTosllBall::initProbMax(){

  EvtId parnum,mesnum,l1num,l2num;
  
  parnum = getParentId();
  mesnum = getDaug(0);
  l1num = getDaug(1);
  l2num = getDaug(2);
  
  //This routine sets the _poleSize.
  double mymaxprob = _calcamp->CalcMaxProb(parnum,mesnum,
					   l1num,l2num,
					   _ballffmodel,_poleSize);

  setProbMax(mymaxprob);

}


void EvtbTosllBall::init(){

  // First choose form factors from the .DEC file
  // 1 = Ali-Ball '01 LCSR
  // 2 = Ali-Ball '99 LCSR
  // 3 = Colangelo 3pt QCD
  // 4 = Melikhov Lattice/Quark dispersion
  // 5 = ???
  // 6 = Ball-Zwicky '05 LCSR (mb = 480)
  // 7 = Ball-Zwicky '05 LCSR (mb = 460 - pseudoscalar modes only)

  // The default is Ali '01
  int theFormFactorModel = 1;

  if(getNArg() == 1)
    theFormFactorModel = (int)getArg(0);

  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be X lepton+ lepton-

  checkSpinParent(EvtSpinType::SCALAR);

  EvtSpinType::spintype mesontype=EvtPDL::getSpinType(getDaug(0));

  if ( !(mesontype == EvtSpinType::VECTOR||
	mesontype == EvtSpinType::SCALAR)) {
    report(Severity::Error,"EvtGen") << "EvtbTosllBall generator expected "
                           << " a SCALAR or VECTOR 1st daughter, found:"<<
                           EvtPDL::name(getDaug(0)).c_str()<<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }

  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::DIRAC);

  _ballffmodel = new EvtbTosllBallFF(theFormFactorModel);
  if (mesontype == EvtSpinType::SCALAR){
    _calcamp = new EvtbTosllScalarAmp(-0.313,4.344,-4.669); 
  } else if (mesontype == EvtSpinType::VECTOR){
    _calcamp = new EvtbTosllVectorAmp(-0.313,4.344,-4.669); 
  }

}
















