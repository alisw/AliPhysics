#ifndef ALIGENEPEMV1_H
#define ALIGENEPEMV1_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * Copyright(c) 1997, 1998, 2002, Adrian Alscher and Kai Hencken          *
 * Copyright(c) 2002 Kai Hencken, Yuri Kharlov, Serguei Sadovsky          *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Event generator of single e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 TeV/nucleon.
// Author: Yuri.Kharlov@cern.ch
// 9 October 2002

#include "AliGenMC.h"
class TEpEmGen;

//-------------------------------------------------------------
class AliGenEpEmv1 : public AliGenMC
{
public:
  AliGenEpEmv1();
  
  virtual ~AliGenEpEmv1();
  void Generate();
  void Init();
  void SetDebug(Int_t debug) {fDebug=debug;}
  
 protected:
  AliGenEpEmv1(const AliGenEpEmv1 & gen);
  AliGenEpEmv1 & operator=(const AliGenEpEmv1 & gen);

  Float_t    fMass;    // electron mass
  TEpEmGen * fEpEmGen; // e+e- generator
  Int_t      fDebug;   // debug level
  Int_t      fEvent;   // internal event number
ClassDef(AliGenEpEmv1,1) // Generator of single e+e- pair production in PbPb ultra-peripheral collisions
};
#endif
