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
// Module: EvtGen/EvtPVVCPLH.hh
//
// Description:
//
// Modification history:
//
//    RYD       November 5, 1999         Module created
//
//    TDP       October 10, 2006         Modified
//
//------------------------------------------------------------------------

#ifndef EVTPVVCPLH_HH
#define EVTPVVCPLH_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtPVVCPLH:public  EvtDecayAmp  {

public:

  EvtPVVCPLH() {}
  virtual ~EvtPVVCPLH();

  virtual std::string getName();
  EvtDecayBase* clone();

  void initProbMax();
  void init();

  void decay(EvtParticle *p); 
  
private:
  bool isBsMixed ( EvtParticle * p ) ;


};

#endif
