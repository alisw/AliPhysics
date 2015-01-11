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
// Module: EvtGen/EvtVVpipi.hh
//
// Description: For decays of a vector to a vector and 2 pions,
//              the decay is assumed to be dominated by S-wave.
//
// Modification history:
//
//    Jim Hunt     June 4, 2008  Module Created
//
//------------------------------------------------------------------------

#ifndef EVTVVPIPI_WEIGHTED_HH
#define EVTVVPIPI_WEIGHTED_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include <string>

class EvtParticle;

class EvtVVPIPI_WEIGHTED:public  EvtDecayAmp  {

public:

  EvtVVPIPI_WEIGHTED() {}
  virtual ~EvtVVPIPI_WEIGHTED();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();

};

#endif
