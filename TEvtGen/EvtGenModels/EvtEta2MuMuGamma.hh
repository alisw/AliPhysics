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
// Module: EvtGen/EvtPi0Dalitz.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTETA2MUMUGAMMA_HH
#define EVTETA2MUMUGAMMA_HH

#include "EvtGenBase/EvtDecayProb.hh"

class EvtParticle;

class EvtEta2MuMuGamma:public  EvtDecayProb  {

public:

  EvtEta2MuMuGamma() {}
  virtual ~EvtEta2MuMuGamma();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();

  void decay(EvtParticle *p);
};

#endif
