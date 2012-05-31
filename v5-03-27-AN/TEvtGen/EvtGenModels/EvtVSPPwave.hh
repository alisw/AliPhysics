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
// Module: EvtGen/EvtVSPPwave.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVSPPWAVE_HH
#define EVTVSPPWAVE_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtVSPPwave:public  EvtDecayAmp  {

public:

  EvtVSPPwave() {}
  virtual ~EvtVSPPwave();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();

};

#endif


