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
// Module: EvtGen/EvtKKLambdaC.hh
//
// Description:Semileptonic decays with pole form form factors
//
// Modification history:
//
//    DJL     April 23, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTKKLAMBDAC_HH
#define EVTKKLAMBDAC_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class Evtparticle;

class EvtKKLambdaC:public  EvtDecayAmp  {

public:

  EvtKKLambdaC() {}
  virtual ~EvtKKLambdaC();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void initProbMax();
  void init();

private:
  EvtSemiLeptonicFF *_ffmodel;
  EvtSemiLeptonicAmp *_calcamp;
};

#endif

