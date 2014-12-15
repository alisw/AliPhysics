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
// Module: EvtGen/EvtSVP.hh
//
// Description:Implementation of the Melikhov semileptonic model
//
// Modification history:
//
//    DJL     April 20, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EvtSVP_HH
#define EvtSVP_HH

#include <fstream>
#include <stdio.h>


#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtParticle;

class EvtSVP:public  EvtDecayAmp  {

public:

  EvtSVP() {}
  virtual ~EvtSVP();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void init();

  virtual void initProbMax();


private:
  int ncall;
};

#endif

