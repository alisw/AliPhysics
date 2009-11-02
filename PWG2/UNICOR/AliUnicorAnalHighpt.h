#ifndef ALIUNICORANALHIGHPT_H
#define ALIUNICORANALHIGHPT_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2008

//=============================================================================
// high-pt correlation analyzer
//=============================================================================

#include "AliUnicorAnal.h"
class AliUnicorEvent;

//=============================================================================
class AliUnicorAnalHighpt : public AliUnicorAnal {
   
 public:
  AliUnicorAnalHighpt(Char_t *nam="highpt", 
	      Double_t emi=-1, Double_t ema=1, 
	      Int_t pid0=0, Int_t pid1=0);   // constructor
  virtual ~AliUnicorAnalHighpt(){}                   // destructor
  void Process(const AliUnicorEvent * const ev0, const AliUnicorEvent * const ev1);  // process 1 (tru) or 2 (mix) events

 protected:
  Int_t    fPid0;                            // particle species 0
  Int_t    fPid1;                            // particle species 1

  ClassDef(AliUnicorAnalHighpt,1)
};
//=============================================================================
#endif
