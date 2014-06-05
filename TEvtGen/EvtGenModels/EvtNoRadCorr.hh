//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2012      University of Warwick, UK
//
// Module: EvtNoRadCorr
//
// Description: Create an empty radiative correction engine which does nothing.
// This is required since the EvtGen constructor still needs at least one 
// concrete implementation of EvtAbsRadCorr for the case when Photos is not used.
//
// Modification history:
//
//    John Back       Sept 2012           Module created
//
//------------------------------------------------------------------------------
//

#ifndef EVTNORADCORR_HH
#define EVTNORADCORR_HH

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include <string>

class EvtParticle;

class EvtNoRadCorr : public EvtAbsRadCorr {

public:
    
  EvtNoRadCorr() {;}
  virtual ~EvtNoRadCorr() {;}

  virtual void doRadCorr(EvtParticle *p) {;}

private:

};

#endif

