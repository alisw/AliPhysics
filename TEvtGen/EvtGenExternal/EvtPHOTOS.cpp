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
// Module: EvtPHOTOS.cc
//
// Description: This routine takes the particle *p and applies
//              the PHOTOS package to generate final state radiation
//              on the produced mesons.
//
// Modification history:
//
//    RYD     October 1, 1997        Module created
//    JJB     May 2011               Modified to use new PHOTOS generator
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenExternal/EvtPHOTOS.hh"
#include "EvtGenExternal/EvtExternalGenFactory.hh"

EvtPHOTOS::EvtPHOTOS() {

  _photosEngine = 0;
  
}

EvtPHOTOS::~EvtPHOTOS() {

}

void EvtPHOTOS::doRadCorr(EvtParticle *p) {

  if (_photosEngine == 0) {
    _photosEngine = EvtExternalGenFactory::getInstance()->getGenerator(EvtExternalGenFactory::PhotosGenId);
  }

  if (_photosEngine != 0) {
    _photosEngine->doDecay(p);
  }
  
}

