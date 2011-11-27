#ifndef ALIPHOSRAWFITTERV1_H
#define ALIPHOSRAWFITTERV1_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id: $ */

// This class extracts the PHOS "digits" of current event
// (amplitude,time, position,gain) from the raw stream 
// provided by AliRawReader. See cxx source for use case.

#include "AliPHOSRawFitterv0.h"
#include "TArrayD.h"
class TList;

class AliPHOSRawFitterv1 : public AliPHOSRawFitterv0 {

public:

  AliPHOSRawFitterv1();
  AliPHOSRawFitterv1(const AliPHOSRawFitterv1& rawFitter);
  AliPHOSRawFitterv1& operator = (const AliPHOSRawFitterv1& rawFitter);
  virtual ~AliPHOSRawFitterv1();
  
  virtual Bool_t Eval(const UShort_t *signal, Int_t sigStart, Int_t sigLength);
  
  static Double_t Gamma2(Double_t dt,Double_t en,Double_t b,TArrayD * fitparams) ; // Shape of correct sample
  //class member function (not object member function)
  void SetLowGainParams(Int_t n, Double_t * params){fSampleParamsLow->Set(n,params) ;}  //fixed parameters of fit function
  void SetHighGainParams(Int_t n,Double_t * params){fSampleParamsHigh->Set(n,params) ;} //fixed parameters of fit function


protected:   
  static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passed to MINUIT


  
protected:  
  TArrayD *fSampleParamsLow ;   //Fixed params of sample parameterization for Low gain 
  TArrayD *fSampleParamsHigh;   //Fixed params of sample parameterization for High gain
  TList   *fToFit ;             //! container to transfer parameters and data to fit



private:
  ClassDef(AliPHOSRawFitterv1,1)
};

#endif
