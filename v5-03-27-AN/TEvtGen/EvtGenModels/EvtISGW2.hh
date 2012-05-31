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
// Module: EvtGen/EvtISGW2.hh
//
// Description:Implementation of the ISGW2 model
// Class to handle semileptonic decays using the ISGW2
// model, as described in PRD 52 5 (1995) by 
// Isgur and Scora.  Electron, muon, and tau models
// are available.  Form factors, q2 and lepton energy
// spectra checked against code from Scora.
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTISGW2_HH
#define EVTISGW2_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtParticle;

class EvtISGW2:public  EvtDecayAmp  {

public:

  EvtISGW2();
  virtual ~EvtISGW2();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void initProbMax();
  void init();

private:
  EvtSemiLeptonicFF *isgw2ffmodel;
  EvtSemiLeptonicAmp *calcamp;
};

#endif

