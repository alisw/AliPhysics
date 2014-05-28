#ifdef EVTGEN_PHOTOS
//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2011      University of Warwick, UK
//
// Module: EvtPhotosEngine
//
// Description: Interface to the PHOTOS external generator
//
// Modification history:
//
//    John Back       May 2011            Module created
//
//------------------------------------------------------------------------

#ifndef EVTPHOTOSENGINE_HH
#define EVTPHOTOSENGINE_HH

#include "EvtGenModels/EvtAbsExternalGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include <string>

class EvtPhotosEngine : public EvtAbsExternalGen {

public:

  EvtPhotosEngine(std::string photonType = "gamma", bool useEvtGenRandom = true);
  virtual ~EvtPhotosEngine();

  virtual bool doDecay(EvtParticle* theMother);

  virtual void initialise();

protected:

private:

  std::string _photonType;
  EvtId _gammaId;
  double _mPhoton;
  bool _initialised;

  HepMC::GenParticle* createGenParticle(EvtParticle* theParticle, bool incoming);

};

#endif

#endif
