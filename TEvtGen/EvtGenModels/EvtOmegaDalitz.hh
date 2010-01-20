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
// Module: EvtGen/EvtOmegaDalitz.hh
//
// Description:Class to handle the omega -> pi pi pi dalitz decay. 

//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTOMEGADALITZ_HH
#define EVTOMEGADALITZ_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtOmegaDalitz:public  EvtDecayAmp  {

public:

  EvtOmegaDalitz() {}
  virtual ~EvtOmegaDalitz();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void decay(EvtParticle *p); 
  void initProbMax();

};

#endif
