#ifndef ITSETFSDD_H
#define ITSETFSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

static const Int_t kMaxNofPoles = 5;
static const Int_t kMaxNofSamples = 1024;

class TString;

class AliITSetfSDD : public TObject {

////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Piergiorgio Cerello
// November 24 1999
//
// AliITSetfSDD is the class describing the electronics for the ITS SDDs. 
//
////////////////////////////////////////////////////////////////////////
  
 public:
    
  AliITSetfSDD() {};                 // default constructor
  AliITSetfSDD(Double_t timestep, Int_t amplif);
  ~AliITSetfSDD() {;}  
  Double_t GetWeightReal(Int_t n) { return fWR[n]; }
  Double_t GetWeightImag(Int_t n) { return fWI[n]; }
  Double_t GetTraFunReal(Int_t n) { return fTfR[n]; }
  Double_t GetTraFunImag(Int_t n) { return fTfI[n]; }
  Int_t GetSamples() { return kMaxNofSamples; }
  Float_t GetTimeDelay() { return fTimeDelay; }
  void PrintElectronics();          // Print Electronics parameters  

 private:

  Float_t  fTimeDelay;         //  Time delay caused by the amplifier shaping
  Double_t fSamplingTime;      //
  Double_t fT0;                //
  Double_t fDf;                //
  Double_t fA0;                //
  Double_t fZeroM[kMaxNofPoles];  // 
  Double_t fZeroR[kMaxNofPoles];  // 
  Double_t fZeroI[kMaxNofPoles];  // 
  Double_t fPoleM[kMaxNofPoles];  // 
  Double_t fPoleR[kMaxNofPoles];  // 
  Double_t fPoleI[kMaxNofPoles];  // 
  Double_t fTfR[kMaxNofSamples];     // Transfer function (real part)
  Double_t fTfI[kMaxNofSamples];     // Transfer function (imaginary part)
  Double_t fWR[kMaxNofSamples];     // Fourier Weights (real part)
  Double_t fWI[kMaxNofSamples];     // Fourier Weights (imaginary part)
  
  ClassDef(AliITSetfSDD,1)  // Class for SDD electornics
    };
    
#endif
  

