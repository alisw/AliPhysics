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
// Module: EvtGen/EvtSVS.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSVS_HH
#define EVTSVS_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;
//Class to handle decays of the form SCALAR->VECTOR SCALAR
//Calls EvtScalarToVectorScalar.

class EvtSVS:public  EvtDecayAmp {

public:

  EvtSVS() {}
  virtual ~EvtSVS();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();

  void initProbMax();

};

#endif
