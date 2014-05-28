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
// Module: EvtGen/EvtTVSPwave.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTTVSPWAVE_HH
#define EVTTVSPWAVE_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

//Class to handles decays of the form TENSOR -> VECTOR SCALAR,
//which proceed via S,P, or D partial waves.  There are
//six arguements, which are the magnetude and phase of
//each partial wave (S, P then D).  Calls EvtTensorToVectorScalar

class EvtTVSPwave:public  EvtDecayAmp  {

public:

  EvtTVSPwave() {}
  virtual ~EvtTVSPwave();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();

  std::string getParamName(int i);
  std::string getParamDefault(int i);
};

#endif
