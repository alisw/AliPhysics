#ifndef ALIMATHBASE_H
#define ALIMATHBASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


 
#include "TObject.h"
#include "TVectorD.h"
#include "TMatrixD.h"

class TH1F;
 
class AliMathBase : public TObject
{
 public:
  AliMathBase();
  virtual ~AliMathBase();
  static void    EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh);
  static void    EvaluateUniExternal(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh, Float_t externalfactor=1);
  static Int_t  Freq(Int_t n, const Int_t *inlist, Int_t *outlist, Bool_t down);    
  static  void TruncatedMean(TH1F * his, TVectorD *param, Float_t down=0, Float_t up=1.0, Bool_t verbose=kFALSE);
  static void LTM(TH1F * his, TVectorD *param=0 , Float_t fraction=1,  Bool_t verbose=kFALSE);
  static Double_t  FitGaus(TH1F* his, TVectorD *param=0, TMatrixD *matrix=0, Float_t xmin=0, Float_t xmax=0,  Bool_t verbose=kFALSE);
    
 ClassDef(AliMathBase,0) // Various mathematical tools for physics analysis - which are not included in ROOT TMath
 
};
#endif
