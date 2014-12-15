//--------------------------------------------------------------------------
//
// Module: EvtGen/EvtHQET2.hh
//
// Description:Implementation of the HQET model with dispersive FF due to 
//             Caprini et al. 
//
// Modification history:
//
//    Marco Bomben   March 11, 2003         Module created
//
//------------------------------------------------------------------------

#ifndef EVTHQET2_HH
#define EVTHQET2_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtParticle;

class EvtHQET2:public  EvtDecayAmp  {

public:

  EvtHQET2();
  virtual ~EvtHQET2();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void initProbMax();
  void init();

private:
  EvtSemiLeptonicFF *hqetffmodel;
  EvtSemiLeptonicAmp *calcamp;
};
#endif



