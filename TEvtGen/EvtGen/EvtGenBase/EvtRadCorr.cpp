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
// Module: EvtRadCorr.cc
//
// Description: RadCorr interface for EvtGen
//              
//
// Modification history:
//
//    Lange     April 27, 2002  - Created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"

#include <stdlib.h>
#include <iostream>
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;


EvtAbsRadCorr* EvtRadCorr::_fsrEngine=0;
bool EvtRadCorr::_alwaysRadCorr=false;
bool EvtRadCorr::_neverRadCorr=false;

EvtRadCorr::EvtRadCorr() {
  _fsrEngine=0;
  _alwaysRadCorr=false;
  _neverRadCorr=false;
}

EvtRadCorr::~EvtRadCorr() {
  if ( _fsrEngine ) delete _fsrEngine;
  _fsrEngine=0;
}

void EvtRadCorr::setRadCorrEngine(EvtAbsRadCorr* fsrEngine){
  _fsrEngine=fsrEngine;
}


void EvtRadCorr::doRadCorr(EvtParticle *p){

  if (_fsrEngine==0){
    report(Severity::Error,"EvtGen") <<"No RadCorr model available in "
			   <<"EvtRadCorr::doRadCorr()."<<endl;
    ::abort();
  }

  if ( !_neverRadCorr) _fsrEngine->doRadCorr(p);
  return;
}


bool EvtRadCorr::alwaysRadCorr() {return _alwaysRadCorr;}
bool EvtRadCorr::neverRadCorr() {return _neverRadCorr;}

void EvtRadCorr::setAlwaysRadCorr() { _alwaysRadCorr=true; _neverRadCorr=false;}
void EvtRadCorr::setNeverRadCorr() { _alwaysRadCorr=false; _neverRadCorr=true;}
void EvtRadCorr::setNormalRadCorr() {_alwaysRadCorr=false; _neverRadCorr=false;}







