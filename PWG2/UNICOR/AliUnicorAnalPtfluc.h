#ifndef ALIUNICORANALPTFLUC_H
#define ALIUNICORANALPTFLUC_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2005

//=============================================================================
// pt-fluctuations analyzer
//=============================================================================

#include "AliUnicorAnal.h"
class AliUnicorEvent;
class AliUnicorHN;

//=============================================================================
class AliUnicorAnalPtfluc : public AliUnicorAnal {
   
 public:
  AliUnicorAnalPtfluc(Char_t *nam="correl", Int_t pid0=0, Int_t pid1=0);    // constructor
  virtual ~AliUnicorAnalPtfluc(){}                                          // destructor
  void Process(Int_t tmr, AliUnicorEvent * const ev0, AliUnicorEvent * const ev1);  // process event(s)

 protected:
  Int_t    fPid0;                       // particle species 0
  Int_t    fPid1;                       // particle species 1
  ClassDef(AliUnicorAnalPtfluc,1)
};
//=============================================================================
#endif
