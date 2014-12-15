#ifndef ALIPHOSRAWFITTERV2_H
#define ALIPHOSRAWFITTERV2_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id: $ */

// This class extracts amplitude, t0 and quality of the PHOS "samples" 
// ising FastFit and two-exponent parameterization

#include "AliPHOSRawFitterv0.h"
class TArrayD ;

class AliPHOSRawFitterv2 : public AliPHOSRawFitterv0 {

public:

  AliPHOSRawFitterv2();
  AliPHOSRawFitterv2(const AliPHOSRawFitterv2& rawFitter);
  AliPHOSRawFitterv2& operator = (const AliPHOSRawFitterv2& rawFitter);
  virtual ~AliPHOSRawFitterv2();

  virtual Bool_t Eval(const UShort_t *signal, Int_t sigStart, Int_t sigLength);
  void SetRawParams(Double_t alpha, Double_t beta){fAlpha=alpha; fBeta=beta;}

private: 
  Bool_t FindAmpT(TArrayD samples, TArrayD times) ;
  void FindMax() ;

private:
  Double_t fAlpha ; //Parameter of sample shape
  Double_t fBeta ;  //Parameter of sample shape
  Double_t fMax ;   //Maximum of parameterization
  
  ClassDef(AliPHOSRawFitterv2,2)
};

#endif
