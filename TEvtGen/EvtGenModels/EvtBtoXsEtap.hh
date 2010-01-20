//--------------------------------------------------------------------------
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
//
// Module: EvtGen/EvtBtoXsEtap.hh
//
// Description:
// Class to generate non-resonant two-body b->s,gluon decays.
//
// Modification history:
//
//    Adlene Hicheur       Jannuary 10, 2001       Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOXSETAP_HH
#define EVTBTOXSETAP_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"

class EvtParticle;

class EvtBtoXsEtap:public  EvtDecayIncoherent  {

public:
  
  EvtBtoXsEtap() {}
  virtual ~EvtBtoXsEtap();

  std::string getName();

  EvtDecayBase* clone();

  void initProbMax();

  void init();

  void decay(EvtParticle *p); 

};

#endif

