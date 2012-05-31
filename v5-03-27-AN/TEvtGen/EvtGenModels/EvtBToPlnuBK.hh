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
// Module: EvtGenModels/EvtBToPlnuBK.hh
//
// Description:   B->Xu l nu with BK (Becirevic-Kaidalov) parametrization
//                Xu is a pseudoscalar (pi_plus,pi0,eta or eta_prime)
//
// Modification history:
//
//    Martin Simard, U. de Montreal, 08/01/2007    Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOPLNUBK_HH
#define EVTBTOPLNUBK_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtParticle;

class EvtBToPlnuBK: public  EvtDecayAmp  {

public:

  EvtBToPlnuBK();
  virtual ~EvtBToPlnuBK();

  std::string getName();
  EvtDecayBase* clone();
  
  void init();
  void initProbMax();

  void decay(EvtParticle *p);

private:

  EvtSemiLeptonicFF *BKmodel;
  EvtSemiLeptonicAmp *calcamp;
  
};

#endif

