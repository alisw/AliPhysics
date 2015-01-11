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
// Module: EvtGen/EvtSVVCP.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSVVCP_HH
#define EVTSVVCP_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtSVVCP:public  EvtDecayAmp  {

public:

  EvtSVVCP() {}
  virtual ~EvtSVVCP();

  std::string getName();
  EvtDecayBase* clone();

  void initProbMax();
  void init();

  void decay(EvtParticle *p); 

  std::string getParamName(int i);
  std::string getParamDefault(int i);
};

#endif
