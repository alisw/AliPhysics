#ifndef ALIPHOSRAWFITTERV4_H
#define ALIPHOSRAWFITTERV4_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id: $ */

// This class extracts the signal parameters (energy, time, quality)
// from ALTRO samples. Energy is in ADC counts, time is in time bin units.
// A coarse algorithm is applied.

class TArrayI;
class TArrayD;
class AliPHOSCalibData;
#include "AliPHOSRawFitterv1.h"

class AliPHOSRawFitterv4 : public AliPHOSRawFitterv1 
{

public:

  AliPHOSRawFitterv4();
  AliPHOSRawFitterv4(const AliPHOSRawFitterv4& rawFitterv4);
  virtual ~AliPHOSRawFitterv4();

  virtual Bool_t Eval(const UShort_t *signal, Int_t sigStart, Int_t sigLength);
  
  //Switch on/off fitting of HighGain samples. By default off
  void FitHighGain(Bool_t on=kTRUE){fFitHighGain=on;}

protected:   
  static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passed to MINUIT
  Bool_t EvalWithFitting(TArrayI*samples, TArrayI * times);

  Bool_t fFitHighGain ; //Switch on/off fitting of the HG channel
  
  ClassDef(AliPHOSRawFitterv4,1)
};

#endif
