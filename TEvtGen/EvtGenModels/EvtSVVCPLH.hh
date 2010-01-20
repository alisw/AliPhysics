//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1999      Caltech, UCSB
//
// Module: EvtGen/EvtSVVCPLH.hh
//
// Description:
//
// Modification history:
//
//    RYD       November 5, 1999         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSVVCPLH_HH
#define EVTSVVCPLH_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtSVVCPLH:public  EvtDecayAmp  {

public:

  EvtSVVCPLH() {}
  virtual ~EvtSVVCPLH();

  std::string getName();
  EvtDecayBase* clone();

  void initProbMax();
  void init();

  void decay(EvtParticle *p); 

};

#endif
