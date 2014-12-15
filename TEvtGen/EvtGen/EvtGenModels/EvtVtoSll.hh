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
// Module: EvtGen/EvtVtoSll.hh
//
// Description:
//
// Modification history:
//
//    RYD    Feb. 28, 2009         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVTOSLL_HH
#define EVTVTOSLL_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtVtoSll:public  EvtDecayAmp  {

public:

  EvtVtoSll() {}
  virtual ~EvtVtoSll();

  std::string getName();
  EvtDecayBase* clone();

  void initProbMax();
  void init();
  void decay(EvtParticle *p); 

};

#endif
