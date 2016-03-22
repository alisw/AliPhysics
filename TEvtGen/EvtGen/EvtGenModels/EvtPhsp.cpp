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
// Module: EvtPhsp.cc
//
// Description: Routine to decay a particle according th phase space
//
// Modification history:
//
//    RYD       January 8, 1997       Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtPhsp.hh"
#include <string>

EvtPhsp::~EvtPhsp() {}

std::string EvtPhsp::getName(){

  return "PHSP";     

}

EvtDecayBase* EvtPhsp::clone(){

  return new EvtPhsp;

}


void EvtPhsp::init(){

  // check that there are 0 arguments
  checkNArg(0);

}

void EvtPhsp::initProbMax(){

  noProbMax();

}

void EvtPhsp::decay( EvtParticle *p ){

  //unneeded - lange - may13-02
  //if ( p->getNDaug() != 0 ) {
    //Will end up here because maxrate multiplies by 1.2
  //  report(Severity::Debug,"EvtGen") << "In EvtPhsp: has "
  //			   <<" daugthers should not be here!"<<endl;
  //  return;
  //}

  p->initializePhaseSpace(getNDaug(),getDaugs());

  return ;
}


