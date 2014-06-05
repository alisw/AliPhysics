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
// Module: EvtGen/EvtPHOTOS.hh
//
// Description: EvtGen's interface to PHOTOS for generation of 
//              QED final state radiation.
//
// Modification history:
//
//    RYD     March 24, 1998         Module created
//    Lange   April 25, 2002 - changed to derive from EvtAbsRadCorr
//
//------------------------------------------------------------------------

#ifndef EVTPHOTOS_HH
#define EVTPHOTOS_HH

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include <string>

class EvtParticle;
class EvtAbsExternalGen;

class EvtPHOTOS : public EvtAbsRadCorr {

public:
    
  EvtPHOTOS();
  virtual ~EvtPHOTOS();

  virtual void doRadCorr(EvtParticle *p);

private:

  EvtAbsExternalGen* _photosEngine;

};

#endif

