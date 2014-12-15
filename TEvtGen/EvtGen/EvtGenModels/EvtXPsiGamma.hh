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
// Module: EvtGen/EvtXPsiGamma.hh
//
// Description:Implementation of the X3872(2-+) -> J/psi gamma decay
//
// Modification history:
//
//    7 May 2012: Module created
//
//------------------------------------------------------------------------

#ifndef EVTXPSIGAMMA_HH
#define EVTXPSIGAMMA_HH

#include <fstream>
#include <stdio.h>

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtParticle;

class EvtXPsiGamma: public EvtDecayAmp {

public:

  EvtXPsiGamma() {}
  virtual ~EvtXPsiGamma();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void init();

  virtual void initProbMax();


private:
//  EvtSemiLeptonicFF *ffmodel;
//  EvtSemiLeptonicAmp *calcamp;
//  int whichfit;
  EvtComplex fT2(EvtVector4R p, EvtVector4R q , EvtTensor4C epsPI, EvtVector4C epsEps, EvtVector4C epsEta); 
  EvtComplex fT3(EvtVector4R p, EvtVector4R q , EvtTensor4C epsPI, EvtVector4C epsEps, EvtVector4C epsEta);
  EvtId _ID0;
  int ncall;
};

#endif

