//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2003      Caltech, UCSB
//
// Module: EvtbTosllAli.cc
//
// Description: Routine to implement b->sll decays according to Ali '02 et al. 
//
// Modification history:
//
//    Ryd     March 30, 2003        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtbTosllAli.hh"
#include "EvtGenModels/EvtbTosllAliFF.hh"
#include "EvtGenModels/EvtbTosllAmp.hh"
#include "EvtGenModels/EvtbTosllScalarAmp.hh"
#include "EvtGenModels/EvtbTosllVectorAmp.hh"

#include <string>
using std::endl;

EvtbTosllAli::~EvtbTosllAli() {}

std::string EvtbTosllAli::getName(){

  return "BTOSLLALI";     
}


EvtDecayBase* EvtbTosllAli::clone(){

  return new EvtbTosllAli;

}

void EvtbTosllAli::decay( EvtParticle *p ){

  setWeight(p->initializePhaseSpace(getNDaug(),getDaugs(),false,
                                    _poleSize,1,2));

  _calcamp->CalcAmp(p,_amp2,_aliffmodel);
  
}


void EvtbTosllAli::initProbMax(){

  EvtId parnum,mesnum,l1num,l2num;
  
  parnum = getParentId();
  mesnum = getDaug(0);
  l1num = getDaug(1);
  l2num = getDaug(2);
  
  //This routine sets the _poleSize.
  double mymaxprob = _calcamp->CalcMaxProb(parnum,mesnum,
					   l1num,l2num,
					   _aliffmodel,_poleSize);

  setProbMax(mymaxprob);

}


void EvtbTosllAli::init(){

  checkNArg(0);
  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be X lepton+ lepton-

  checkSpinParent(EvtSpinType::SCALAR);

  EvtSpinType::spintype mesontype=EvtPDL::getSpinType(getDaug(0));

  if ( !(mesontype == EvtSpinType::VECTOR||
	mesontype == EvtSpinType::SCALAR)) {
    report(Severity::Error,"EvtGen") << "EvtbTosllAli generator expected "
                           << " a SCALAR or VECTOR 1st daughter, found:"<<
                           EvtPDL::name(getDaug(0)).c_str()<<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }

  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::DIRAC);

  _aliffmodel = new EvtbTosllAliFF();
  if (mesontype == EvtSpinType::SCALAR){
    _calcamp = new EvtbTosllScalarAmp(-0.313,4.344,-4.669); 
  }
  if (mesontype == EvtSpinType::VECTOR){
    _calcamp = new EvtbTosllVectorAmp(-0.313,4.344,-4.669); 
  }

}





