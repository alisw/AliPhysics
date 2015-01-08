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
// Module: EvtGen/EvtVVSPwave.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVVSPWAVE_HH
#define EVTVVSPWAVE_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtVVSPwave:public  EvtDecayAmp  {

public:

  EvtVVSPwave() {}
  virtual ~EvtVVSPwave();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();

};

#endif
