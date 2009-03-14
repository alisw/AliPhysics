#ifndef ALITPCKALMANDEDX_H
#define ALITPCKALMANDEDX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TNamed.h"
#include "TMatrixD.h"
#include "TVectorD.h"
class TTreeSRedirector;
class TObjArray;

class AliTPCkalmandEdx: public TNamed{
public:
  AliTPCkalmandEdx();
  AliTPCkalmandEdx(const char* name, const char* title, Int_t sampleSize=50);
  AliTPCkalmandEdx(const AliTPCkalmandEdx & kalman);
  void UpdatedEdxPair(Int_t ip0, Int_t ip1,
		      Double_t dedx0, Double_t dedx1, 
		      Double_t s0, Double_t s1, 
		      Double_t kY0, Double_t kY1,
		      Double_t kZ0, Double_t kZ1,
		      Double_t dR0, Double_t dR1, 
		      TTreeSRedirector *debug=0);
  void UpdatedEdx(Int_t ip0,
		  Double_t dedx0, 
		  Double_t dedxRef, 
		  Double_t s0,
		  Double_t kY0,
		  Double_t kZ0,
		  Double_t dR0,
		  TTreeSRedirector *debug=0);
  static Int_t GetIndex(Int_t i, Int_t j) { return 15*i+j;}
  //
public:
  void AdddEdx(Int_t ip0,Double_t dedx0, Double_t dedxRef);
  void Init();
  TMatrixD * fState;           // state vector
  TMatrixD * fCovariance;      // covariance
  TMatrixD * fMat1;            // helper unit matrix
  Int_t      fNpad;            // number of pad types
  Int_t      fNpar;            // number of parameters
  Int_t      fNelem;           // number of elements
  //
  // initial parameters estimate
  //
  Int_t      fSampleSize;       // size of  starting  sample
  Int_t      fInit;             // number of initialized estimators
  //
  TVectorD   fSample[3];        // !training sample for robust estimate of initial parameters
  TVectorD   fSampleStat[3];    // sample statistic
  Int_t      fCounter[3];       // counter of samples
  //
  ClassDef(AliTPCkalmandEdx,1);
};

#endif

