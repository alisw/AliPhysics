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
// Module: EvtGen/EvtISGW.hh
//
// Description:Implementation of the ISGW model
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTISGW_HH
#define EVTISGW_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"
class EvtParticle;



class EvtISGW:public  EvtDecayAmp  {

public:

  EvtISGW() {}
  virtual ~EvtISGW();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void init();

private:
  EvtSemiLeptonicFF *isgwffmodel;
  EvtSemiLeptonicAmp *calcamp;
};



#endif

