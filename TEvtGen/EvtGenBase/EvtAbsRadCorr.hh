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
// Module: EvtGenBase/EvtAbsRadCorr.hh
//
// Description: 
//
// Modification history:
//
//    Lange April 25, 2002   - Created
//
//------------------------------------------------------------------------

#ifndef EVTABSRADCORR_HH
#define EVTABSRADCORR_HH

#include <assert.h>
#include <iostream>
class EvtParticle;

class EvtAbsRadCorr {

public:

  EvtAbsRadCorr() {};
  virtual ~EvtAbsRadCorr() {};
  virtual void doRadCorr(EvtParticle *p)=0;


private:

};

#endif
