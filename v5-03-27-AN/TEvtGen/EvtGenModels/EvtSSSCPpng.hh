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
// Module: EvtGen/EvtSSSCPpng.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSSSCPPNG_HH
#define EVTSSSCPPNG_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtSSSCPpng:public  EvtDecayAmp  {

public:

  EvtSSSCPpng() {}
  virtual ~EvtSSSCPpng();

  std::string getName();
  EvtDecayBase* clone();

  void initProbMax(); 

  void init();
  void decay(EvtParticle *p); 

};

#endif
