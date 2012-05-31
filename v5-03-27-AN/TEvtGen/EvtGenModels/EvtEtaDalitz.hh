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
// Module: EvtGen/EvtEtaDalitz.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTETADALITZ_HH
#define EVTETADALITZ_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtEtaDalitz:public  EvtDecayAmp  {

public:

  EvtEtaDalitz() {}
  virtual ~EvtEtaDalitz();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();

  void decay(EvtParticle *p); 

};

#endif
