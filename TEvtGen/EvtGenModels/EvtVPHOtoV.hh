//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2004      Cornell
//
// Module: EvtGen/EvtVPHOtoV.hh
//
// Description:
//
// Modification history:
//
//    Ryd      March 9, 2004      Module created
//
//------------------------------------------------------------------------

#ifndef EVTVPHOTOV_HH
#define EVTVPHOTOV_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtVPHOtoV:public  EvtDecayAmp  {

public:

  EvtVPHOtoV() {}
  virtual ~EvtVPHOtoV();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();

};

#endif
