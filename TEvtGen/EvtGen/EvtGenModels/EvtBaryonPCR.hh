//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtBaryonPCR.hh
//
// Description:Implementation of the BaryonPCR model
// Class to handle semileptonic decays using the BaryonVminusA
// model.
//
// Modification history:
//
//    R.J. Tesarek     May 28, 2004     Module created
//    Karen Gibson     1/20/2006        Module updated for 1/2+->1/2+,
//                                      1/2+->1/2-, 1/2+->3/2- Lambda decays
//
//------------------------------------------------------------------------

#ifndef EVTBARYONPCR_HH
#define EVTBARYONPCR_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicBaryonAmp.hh"

class EvtParticle;

class EvtBaryonPCR:public  EvtDecayAmp  {

public:

  EvtBaryonPCR();
  virtual ~EvtBaryonPCR();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void initProbMax();
  void init();

private:
  EvtSemiLeptonicFF *baryonpcrffmodel;
  EvtSemiLeptonicBaryonAmp *calcamp;
};

#endif

