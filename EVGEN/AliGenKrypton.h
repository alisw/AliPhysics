#ifndef ALIGENKR_H
#define ALIGENKR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Generator of Krypton decay
//
#include "AliGenerator.h"
class AliGenKrypton : public AliGenerator
{
public:
  AliGenKrypton();
  virtual void Generate();
  virtual ~AliGenKrypton(){}
  private:
  void KrDecay(Int_t &nelectron, Int_t &ngamma, Double_t *eelectron, Double_t *egamma);

  ClassDef(AliGenKrypton,1)
};
#endif
