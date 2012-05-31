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
// Module: EvtSLBKPole.cc
//
// Description: Routine to implement semileptonic decays according
//              to light cone sum rules
//
// Modification history:
//
//    liheng       October 20, 2005       Module created
//
//------------------------------------------------------------------------
// 
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtSLBKPole.hh"//modified
#include "EvtGenModels/EvtSLBKPoleFF.hh"//modified
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicTensorAmp.hh"
#include <string>

EvtSLBKPole::~EvtSLBKPole() {}

std::string EvtSLBKPole::getName(){

   return "SLBKPOLE";//modified

}


EvtDecayBase* EvtSLBKPole::clone(){//modified

  return new EvtSLBKPole;

}

void EvtSLBKPole::decay( EvtParticle *p ){//modified

  p->initializePhaseSpace(getNDaug(),getDaugs());

  calcamp->CalcAmp(p,_amp2,SLBKPoleffmodel);//modified
  return;
}

void EvtSLBKPole::initProbMax(){

EvtId parnum,mesnum,lnum,nunum;

parnum = getParentId();
mesnum = getDaug(0);
lnum = getDaug(1);
nunum = getDaug(2);

double mymaxprob = calcamp->CalcMaxProb(parnum,mesnum,
			lnum,nunum,SLBKPoleffmodel);//modified

setProbMax(mymaxprob);

}


void EvtSLBKPole::init(){//modified
  
  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be X lepton neutrino

  checkSpinParent(EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);

  EvtSpinType::spintype mesontype=EvtPDL::getSpinType(getDaug(0));

  SLBKPoleffmodel = new EvtSLBKPoleFF(getNArg(),getArgs());//modified
  
  if ( mesontype==EvtSpinType::SCALAR ) { 
    calcamp = new EvtSemiLeptonicScalarAmp; 
  }
  if ( mesontype==EvtSpinType::VECTOR ) { 
    calcamp = new EvtSemiLeptonicVectorAmp; 
  }
  if ( mesontype==EvtSpinType::TENSOR ) { 
    calcamp = new EvtSemiLeptonicTensorAmp; 
  }
  
}

