#ifndef ALIPHOSRAWFITTERV2_H
#define ALIPHOSRAWFITTERV2_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id: $ */

// This class extracts the PHOS "digits" of current event
// (amplitude,time, position,gain) from the raw stream 
// provided by AliRawReader. See cxx source for use case.

#include "AliPHOSRawFitterv0.h"
class TList;

class AliPHOSRawFitterv2 : public AliPHOSRawFitterv0 {

public:

  AliPHOSRawFitterv2();
  AliPHOSRawFitterv2(const AliPHOSRawFitterv2& rawFitter);
  AliPHOSRawFitterv2& operator = (const AliPHOSRawFitterv2& rawFitter);
  virtual ~AliPHOSRawFitterv2();

  virtual Bool_t Eval(const UShort_t *signal, Int_t sigStart, Int_t sigLength);

  void SetNTimeSamples(Short_t n=25)     { fNtimeSamples=n ;}
  void SetLowGainTParams (Double_t *pars){ for(Int_t i=0;i<3;i++) fLGpar[i]=pars[i] ;}
  void SetHighGainTParams(Double_t *pars){ for(Int_t i=0;i<3;i++) fHGpar[i]=pars[i] ;}
  void SetRMScut(Double_t cut=2.)        { fRMScut = cut ;}

private:
  Short_t  fNtimeSamples ; //Number of samples (after start) used to extract time
  Double_t fLGpar[3] ;     //parameters for shape parameterization
  Double_t fHGpar[3] ;     //parameters for shape parameterization
  Double_t fRMScut ;       //cut to estmate goodness of sample
  
  ClassDef(AliPHOSRawFitterv2,1)
};

#endif
