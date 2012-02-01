#ifndef ALIUNICORANALSINGLE_H
#define ALIUNICORANALSINGLE_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// single particle analyzer
//=============================================================================

#include "AliUnicorAnal.h"
class AliUnicorEvent;
class AliUnicorHN;

//=============================================================================
class AliUnicorAnalSingle : public AliUnicorAnal {
   
 public:
  AliUnicorAnalSingle(const char *nam="single", 
	      Double_t emi=-1, Double_t ema=1, 
	      Int_t pid=0);                       // constructor
  virtual ~AliUnicorAnalSingle(){}                        // destructor
  void Process(AliUnicorEvent *ev);                       // fill histograms

 protected:
  Int_t    fPid;                                  // pid; 0 means all
  Double_t fMass;                                 // mass (if pid!=0)

  ClassDef(AliUnicorAnalSingle,1)
};
//=============================================================================
#endif
