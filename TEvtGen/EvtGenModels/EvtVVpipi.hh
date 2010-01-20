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
//    RYD/CAHN  December 10, 1999         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVVPIPI_HH
#define EVTVVPIPI_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtVVpipi:public  EvtDecayAmp  {

public:

  EvtVVpipi() {}
  virtual ~EvtVVpipi();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();

};

#endif

