//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2002      Caltech
//
// Module: EvtGen/EvtKstarstargamma.hh
//
// Description:
//
// Modification history:
//
//    RYD     Nov 15, 2002         Module created
//
//------------------------------------------------------------------------

#ifndef EVTKSTARSTARGAMMA_HH
#define EVTKSTARSTARGAMMA_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtKstarstargamma:public  EvtDecayAmp  {

public:

  EvtKstarstargamma() {}
  virtual ~EvtKstarstargamma();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();

};

#endif
