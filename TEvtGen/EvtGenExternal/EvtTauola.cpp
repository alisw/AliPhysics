//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C)  2011      University of Warwick, UK
//
// Description: Use the Tauola external generator for tau decays
//
// Modification history:
//
//    John Back    May 2011   Module created
//
//------------------------------------------------------------------------

#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenExternal/EvtTauola.hh"

#include "EvtGenExternal/EvtExternalGenFactory.hh"
#include "EvtGenModels/EvtAbsExternalGen.hh"
#include "EvtGenBase/EvtDecayBase.hh"

#include <string>
#include <iostream>
#include <cmath>

EvtTauola::EvtTauola() {

  // Set the Tauola engine to a null pointer at first.
  // When we do the decay, we retrieve the pointer to the Tauola engine
  // and use that for all decays. All clones will use the same Tauola engine.
  _tauolaEngine = 0;

}

EvtTauola::~EvtTauola() {}

std::string EvtTauola::getName(){

  return "TAUOLA";

}

EvtDecayBase* EvtTauola::clone(){

  return new EvtTauola();

}


void EvtTauola::init(){


}

void EvtTauola::initProbMax(){

  noProbMax();

}

void EvtTauola::decay( EvtParticle *p ){

  // We check to see if the Tauola engine has been created before doing the decay.
  // This should only create the full Tauola engine once, and all clones will 
  // point to the same engine.

  if (_tauolaEngine == 0) {
    _tauolaEngine = EvtExternalGenFactory::getInstance()->getGenerator(EvtExternalGenFactory::TauolaGenId);
  }
    
  if (_tauolaEngine != 0) {
    _tauolaEngine->doDecay(p);
  }

}
