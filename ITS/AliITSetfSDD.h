#ifndef ITSETFSDD_H
#define ITSETFSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class AliITSetfSDD : public TObject {

////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Piergiorgio Cerello
// November 24 1999
//
// AliITSetfSDD is the class describing the electronics for the ITS SDDs. 
//
// Data members:
//
////////////////////////////////////////////////////////////////////////

 private:
  
  Double_t fSamplingTime;      //
  Double_t fT0;                //
  Double_t fDf;                //
  Double_t fA0;                //
  static const Int_t fMaxNofPoles = 5;
  Double_t fZero_M[fMaxNofPoles];  // 
  Double_t fZero_R[fMaxNofPoles];  // 
  Double_t fZero_I[fMaxNofPoles];  // 
  Double_t fPole_M[fMaxNofPoles];  // 
  Double_t fPole_R[fMaxNofPoles];  // 
  Double_t fPole_I[fMaxNofPoles];  // 
  static const Int_t fMaxNofSamples = 256;
  Double_t fTf_R[fMaxNofSamples];     // Transfer function (real part)
  Double_t fTf_I[fMaxNofSamples];     // Transfer function (imaginary part)
  Double_t fW_R[fMaxNofSamples];     // Fourier Weights (real part)
  Double_t fW_I[fMaxNofSamples];     // Fourier Weights (imaginary part)
  
 public:
    
  AliITSetfSDD() {};                 // default constructor
  AliITSetfSDD(Double_t);
  ~AliITSetfSDD() {;}  
  Double_t GetWeightReal(Int_t n) { return fW_R[n]; }
  Double_t GetWeightImag(Int_t n) { return fW_I[n]; }
  Double_t GetTraFunReal(Int_t n) { return fTf_R[n]; }
  Double_t GetTraFunImag(Int_t n) { return fTf_I[n]; }
  Int_t GetSamples() { return fMaxNofSamples; }
  void Print();          // Print Electronics parameters  
  
  friend class AliITSmapSDD;
  ClassDef(AliITSetfSDD,1)  // Class for SDD electornics
    };
    
#endif
  

