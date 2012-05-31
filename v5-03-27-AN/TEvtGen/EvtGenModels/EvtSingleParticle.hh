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
// Module: EvtGen/EvtSingleParticle.hh
//
// Description:
//This is a special decay model to generate single particles.
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------


#ifndef EVTSINGLEPARTICLE_HH
#define EVTSINGLEPARTICLE_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"
class EvtParticle;

class EvtSingleParticle:public  EvtDecayIncoherent  {

public:

  EvtSingleParticle() {}
  virtual ~EvtSingleParticle();

  std::string getName();
  
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 

  void init();

private:

  double pmin,pmax;
  double cthetamin,cthetamax;
  double phimin,phimax;

};

#endif

