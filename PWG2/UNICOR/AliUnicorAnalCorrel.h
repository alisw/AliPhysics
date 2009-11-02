#ifndef ALIUNICORANALCORREL_H
#define ALIUNICORANALCORREL_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2005

//=============================================================================
// two-particle correlation analyzer
//=============================================================================

#include "AliUnicorAnal.h"
#include "AliUnicorPair.h"
class AliUnicorEvent;

//=============================================================================
class AliUnicorAnalCorrel : public AliUnicorAnal {
   
 public:
  AliUnicorAnalCorrel(Char_t *nam="correl", 
	      Double_t emi=-1, Double_t ema=1, 
	      Int_t pid0=0, Int_t pid1=0);     // constructor
  virtual ~AliUnicorAnalCorrel(){}                     // destructor
  // process one (tru) or two (mix) events
  void Process(Int_t tmr, const AliUnicorEvent * const ev0, const AliUnicorEvent * const ev1, Double_t phirot);

 protected:
  Int_t    fPid0;                       // particle species 0
  Int_t    fPid1;                       // particle species 1
  Double_t fMass0;                      // mass 0
  Double_t fMass1;                      // mass 1
  AliUnicorPair    fPa;                         // pair buffer for calculations

  ClassDef(AliUnicorAnalCorrel,1)
};
//=============================================================================
#endif
