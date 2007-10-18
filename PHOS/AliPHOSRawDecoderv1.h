#ifndef ALIPHOSRAWDECODERV1_H
#define ALIPHOSRAWDECODERV1_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */

// This class extracts the PHOS "digits" of current event
// (amplitude,time, position,gain) from the raw stream 
// provided by AliRawReader. See cxx source for use case.

#include "AliRawReader.h"
#include "AliCaloRawStream.h"
#include "AliPHOSRawDecoder.h"


class AliPHOSRawDecoderv1 : public AliPHOSRawDecoder {

public:

  AliPHOSRawDecoderv1();
  AliPHOSRawDecoderv1(AliRawReader* rawReader, AliAltroMapping **mapping = NULL);
  AliPHOSRawDecoderv1(const AliPHOSRawDecoderv1& rawDecoder);
  AliPHOSRawDecoderv1& operator = (const AliPHOSRawDecoderv1& rawDecoder);
  virtual ~AliPHOSRawDecoderv1();

  virtual Bool_t NextDigit();

  static Double_t Gamma2(Double_t dt,Double_t p,Double_t en,Double_t a) ; // Shape of correct sample
                                                 //class member function (not object member function)
  static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passed to MINUIT
private:
  
  ClassDef(AliPHOSRawDecoderv1,1)
};

#endif
