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
////////////////////////////////////////////////////////////////////////
  
 public:
    
  AliITSetfSDD() {};                 // default constructor
  AliITSetfSDD(Double_t timestep);
  ~AliITSetfSDD() {;}  
  Double_t GetWeightReal(Int_t n) { return fWR[n]; }
  Double_t GetWeightImag(Int_t n) { return fWI[n]; }
  Double_t GetTraFunReal(Int_t n) { return fTfR[n]; }
  Double_t GetTraFunImag(Int_t n) { return fTfI[n]; }
  Int_t GetSamples() { return fkMaxNofSamples; }
  void PrintElectronics();          // Print Electronics parameters  

 private:

  static const Int_t fkMaxNofPoles = 5;
  static const Int_t fkMaxNofSamples = 1024;
  
  Double_t fSamplingTime;      //
  Double_t fT0;                //
  Double_t fDf;                //
  Double_t fA0;                //
  Double_t fZeroM[fkMaxNofPoles];  // 
  Double_t fZeroR[fkMaxNofPoles];  // 
  Double_t fZeroI[fkMaxNofPoles];  // 
  Double_t fPoleM[fkMaxNofPoles];  // 
  Double_t fPoleR[fkMaxNofPoles];  // 
  Double_t fPoleI[fkMaxNofPoles];  // 
  Double_t fTfR[fkMaxNofSamples];     // Transfer function (real part)
  Double_t fTfI[fkMaxNofSamples];     // Transfer function (imaginary part)
  Double_t fWR[fkMaxNofSamples];     // Fourier Weights (real part)
  Double_t fWI[fkMaxNofSamples];     // Fourier Weights (imaginary part)
  
  ClassDef(AliITSetfSDD,1)  // Class for SDD electornics
    };
    
#endif
  

