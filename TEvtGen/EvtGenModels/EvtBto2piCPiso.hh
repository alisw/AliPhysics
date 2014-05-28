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
// Module: EvtGen/EvtBto2piCPiso.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTO2PICPISO_HH
#define EVTBTO2PICPISO_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtParticle.hh"

class EvtBto2piCPiso:public  EvtDecayAmp  {

public:
  
  EvtBto2piCPiso() {}
  virtual ~EvtBto2piCPiso();
  
  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();

  void decay(EvtParticle *p); 

  std::string getParamName(int i);
  std::string getParamDefault(int i);
};

#endif
