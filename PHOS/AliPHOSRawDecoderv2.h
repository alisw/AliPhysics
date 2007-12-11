#ifndef ALIPHOSRAWDECODERV2_H
#define ALIPHOSRAWDECODERV2_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */

// This class extracts the PHOS "digits" of current event
// (amplitude,time, position,gain) from the raw stream 
// provided by AliRawReader. See cxx source for use case.

#include "AliRawReader.h"
#include "AliCaloRawStream.h"
#include "AliPHOSRawDecoder.h"
#include "TArrayD.h"
class TList;

class AliPHOSRawDecoderv2 : public AliPHOSRawDecoder {

public:

  AliPHOSRawDecoderv2();
  AliPHOSRawDecoderv2(AliRawReader* rawReader, AliAltroMapping **mapping = NULL);
  AliPHOSRawDecoderv2(const AliPHOSRawDecoderv2& rawDecoder);
  AliPHOSRawDecoderv2& operator = (const AliPHOSRawDecoderv2& rawDecoder);
  virtual ~AliPHOSRawDecoderv2();

  virtual Bool_t NextDigit();

  void SetNTimeSamples(Short_t n=25){fNtimeSamples=n ;} 
  void SetLowGainTParams(Double_t *pars){ for(Int_t i=0;i<3;i++) fLGpar[i]=pars[i] ; }  
  void SetHighGainTParams(Double_t *pars){ for(Int_t i=0;i<3;i++) fHGpar[i]=pars[i] ; }  
  void SetRMScut(Double_t cut=2.){fRMScut = cut ;}
private:
  Short_t fNtimeSamples ;  //Number of samples (after start) used to extract time
  Double_t fLGpar[3] ;     //parameters for shape parameterization
  Double_t fHGpar[3] ;     //parameters for shape parameterization
  Double_t fRMScut ;       //cut to estmate goodness of sample
  
  ClassDef(AliPHOSRawDecoderv2,1)
};

#endif
