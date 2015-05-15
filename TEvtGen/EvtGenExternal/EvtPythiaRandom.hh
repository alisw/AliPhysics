#ifdef EVTGEN_PYTHIA
//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2013    University of Warwick, UK
//
// Module: EvtGenExternal/EvtPythiaRandom.hh
//
// Description: Class to specify the chosen EvtGen random number (engine)
// to also be used for Pythia 8.
//
// Modification history:
//
//    JJB     January 2013      Module created
//
//------------------------------------------------------------------------

#ifndef EVTPYTHIARANDOM_HH
#define EVTPYTHIARANDOM_HH

#include "EvtGenBase/EvtRandom.hh"

#include "Pythia8/Basics.h"

class EvtPythiaRandom : public Pythia8::RndmEngine {

public:

  EvtPythiaRandom() {};

  virtual ~EvtPythiaRandom() {};

  virtual double flat() {return EvtRandom::Flat();}

private:

};

#endif

#endif
