//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2011      University of Warwick, UK
//
// Module: EvtAbsExternalGen
//
// Description: Pure abstract interface for external physics generators
//
// Modification history:
//
//    John Back       April 2011            Module created
//
//------------------------------------------------------------------------

#ifndef EVTABS_EXTERNALGEN_HH
#define EVTABS_EXTERNALGEN_HH

#include "EvtGenBase/EvtParticle.hh"

class EvtAbsExternalGen {

public:

  EvtAbsExternalGen() {};
  virtual ~EvtAbsExternalGen() {};

  virtual bool doDecay(EvtParticle* theMother) = 0;
  virtual void initialise() = 0;

protected:

private:

};

#endif
