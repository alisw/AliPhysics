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
// Module: EvtGen/EvtBTo3piCP.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTO3PICP_HH
#define EVTBTO3PICP_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;



class EvtBTo3piCP:public  EvtDecayAmp  {

public:

  EvtBTo3piCP() {}
  virtual ~EvtBTo3piCP();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();

  void decay(EvtParticle *p); 

};

#endif
