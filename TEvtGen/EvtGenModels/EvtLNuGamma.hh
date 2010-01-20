//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2001  Caltech
//
// Module: EvtGen/EvtLNuGamma.hh
//
// Description: B+ -> l+ nu gamma.  Form factor is tree level, from 
//  Korchemsky, Pirjol, and Yan,Phy Rev D 61 (200) 114510
//               
//
// Modification history:
//
//    Edward Chen     April 24, 2001         Module created
//
//------------------------------------------------------------------------

#ifndef EVTLNUGAMMA_HH
#define EVTLNUGAMMA_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtLNuGamma:public  EvtDecayAmp  {

public:

  EvtLNuGamma();
  virtual ~EvtLNuGamma();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();
  double getFormFactor(double photonEnergy);

  bool _fafvzero;

};

#endif
