#ifndef ALIMATHBASE_H
#define ALIMATHBASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


 
#include "TObject.h"

 
class AliMathBase : public TObject
{
 public:
  AliMathBase();
  virtual ~AliMathBase();
  static void    EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh);
  static void    EvaluateUniExternal(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh, Float_t externalfactor=1);
  static Int_t  Freq(Int_t n, const Int_t *inlist, Int_t *outlist, Bool_t down);    

 ClassDef(AliMathBase,0) // Various mathematical tools for physics analysis - which are not included in ROOT TMath
 
};
#endif
