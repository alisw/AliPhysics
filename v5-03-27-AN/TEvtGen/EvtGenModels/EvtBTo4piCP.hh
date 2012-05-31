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
// Module: EvtGen/EvtBTo4piCP.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTO4PICP_HH
#define EVTBTO4PICP_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtBTo4piCP:public  EvtDecayAmp  {

public:

  EvtBTo4piCP() {}
  virtual ~EvtBTo4piCP();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void decay(EvtParticle *p); 

};

#endif
