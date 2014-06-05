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
// Module: EvtDecayParm.cc
//
// Description: Store decay parameters for one decay.
//
// Modification history:
//
//    RYD     April 5, 1997         Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>
#include "EvtGenBase/EvtParticleDecay.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include <vector>
using std::fstream;
void EvtParticleDecay::printSummary(){

  if (_decay!=0) _decay->printSummary();

}

void EvtParticleDecay::chargeConj(EvtParticleDecay *decay){

  _brfrsum=decay->_brfrsum;
  _massmin=decay->_massmin;

  _decay=decay->_decay->clone();

  int ndaug=decay->_decay->getNDaug();
  int narg=decay->_decay->getNArg();
  double brfr=decay->_decay->getBranchingFraction();
  std::string name=decay->_decay->getName();
  EvtId ipar=EvtPDL::chargeConj(decay->_decay->getParentId());
  int i;
  EvtId* daug=new EvtId[ndaug];
  for(i=0;i<ndaug;i++){
    daug[i]=EvtPDL::chargeConj(decay->_decay->getDaug(i));
  }
  //Had to add 1 to make sure the vector is not empty!
  std::vector<std::string> args;
  for(i=0;i<narg;i++){
    args.push_back(decay->_decay->getArgStr(i));
  }

  _decay->saveDecayInfo(ipar,ndaug,daug,narg,args,name,brfr);

  if (decay->_decay->getPHOTOS()) _decay->setPHOTOS();

  delete [] daug;

}







