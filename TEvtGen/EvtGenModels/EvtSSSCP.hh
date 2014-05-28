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
// Module: EvtGen/EvtSSSCP.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSSSCP_HH
#define EVTSSSCP_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;


class EvtSSSCP:public  EvtDecayAmp  {

public:

  EvtSSSCP() {}
  virtual ~EvtSSSCP();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();

  void decay(EvtParticle *p); 

  std::string getParamName(int i);
  std::string getParamDefault(int i);
};

#endif
