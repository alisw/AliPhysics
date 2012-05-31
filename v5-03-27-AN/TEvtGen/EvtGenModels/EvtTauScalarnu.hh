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
// Module: EvtGen/EvtTauScalarnu.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTTAUSCALARNU_HH
#define EVTTAUSCALARNU_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtTauScalarnu:public  EvtDecayAmp  {

public:
  
  EvtTauScalarnu() {}
  virtual ~EvtTauScalarnu();

  std::string getName();
  EvtDecayBase* clone();

  void initProbMax();
  void init();
  void decay(EvtParticle *p); 

};

#endif
