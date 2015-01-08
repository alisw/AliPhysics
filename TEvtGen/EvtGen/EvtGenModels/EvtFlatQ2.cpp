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
// Module: EvtFlatQ2.cc
//
// Description: B->Xu l nu with flat q2 distribution
//
// Modification history:
//
//    David Cote, U. de Montreal, 11/02/2003    Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenModels/EvtFlatQ2.hh"

#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <string>
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
using std::fstream;

EvtFlatQ2::~EvtFlatQ2() {}

std::string EvtFlatQ2::getName(){

    return "FLATQ2";

}

EvtDecayBase* EvtFlatQ2::clone(){

  return new EvtFlatQ2;

}


void EvtFlatQ2::initProbMax(){

  setProbMax(100);

}


void EvtFlatQ2::init(){

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(3);

  //We expect B->X l nu events
  checkSpinParent(EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);

}


void EvtFlatQ2::decay( EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtVector4R p4Xu = p->getDaug(0)->getP4();
  double pXu_x2=p4Xu.get(1)*p4Xu.get(1);
  double pXu_y2=p4Xu.get(2)*p4Xu.get(2);
  double pXu_z2=p4Xu.get(3)*p4Xu.get(3);
  double pXu = sqrt(pXu_x2+pXu_y2+pXu_z2);
  double prob=1/pXu;
  
  if(pXu>0.01) setProb(prob);

  return;
}


