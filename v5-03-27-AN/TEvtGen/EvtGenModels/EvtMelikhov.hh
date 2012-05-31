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
// Module: EvtGen/EvtMelikhov.hh
//
// Description:Implementation of the Melikhov semileptonic model
//
// Modification history:
//
//    DJL     April 20, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTMELIKHOV_HH
#define EVTMELIKHOV_HH

#include <fstream>
#include <stdio.h>


#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtParticle;

class EvtMelikhov:public  EvtDecayAmp  {

public:

  EvtMelikhov() {}
  virtual ~EvtMelikhov();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void init();

private:
  EvtSemiLeptonicFF *Melikhovffmodel;
  EvtSemiLeptonicAmp *calcamp;
};

#endif

