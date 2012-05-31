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
// Module: EvtGen/EvtBHadronic.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------
#ifndef EVTBHADRONIC_HH
#define EVTBHADRONIC_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtBHadronic:public  EvtDecayAmp{

public:

  EvtBHadronic() {}
  virtual ~EvtBHadronic();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void decay(EvtParticle *p); 

};

#endif
