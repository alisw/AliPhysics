#ifndef ALIITSETFSDD_H
#define ALIITSETFSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////////////////////////////
// Version: 0                                                         //
// Written by Piergiorgio Cerello                                     //
// November 24 1999                                                   //  
//                                                                    //
// AliITSetfSDD is the class describing                               //
// the electronics for the ITS SDDs.                                  //
////////////////////////////////////////////////////////////////////////

#include <TObject.h>


class TString;

class AliITSetfSDD : public TObject {

  
 public:
    
  AliITSetfSDD();                 // default constructor
  AliITSetfSDD(Double_t timestep, Int_t amplif); 
  virtual ~AliITSetfSDD();  
  Double_t GetWeightReal(Int_t n) const { return fWR[n]; }
  Double_t GetWeightImag(Int_t n) const { return fWI[n]; }
  Double_t GetTraFunReal(Int_t n) const { return fTfR[n]; }
  Double_t GetTraFunImag(Int_t n) const { return fTfI[n]; }
  Float_t GetTimeDelay() const { return fTimeDelay; }
  void PrintElectronics() const; // Print Electronics parameters  
  class AliITSetfSDDparam {
    public :
      static Int_t NumberOfPoles() {return fgkMaxNofPoles;}
      static Int_t NumberOfSamples(){return fgkMaxNofSamples;}
  private:
      static const Int_t fgkMaxNofPoles; // Set to 5  in impl. file
      static const Int_t fgkMaxNofSamples; // Set to 1024 in impl. file
  };
  static Int_t GetSamples() { return AliITSetfSDDparam::NumberOfSamples(); }
 protected:

  // copy constructor (NO copy allowed: the constructor is protected
  // to avoid misuse)
  AliITSetfSDD(const AliITSetfSDD& obj);
  // assignment operator (NO assignment allowed)
  AliITSetfSDD& operator=(const AliITSetfSDD& /* obj */);

  Float_t  fTimeDelay;         //  Time delay caused by the amplifier shaping
  Double_t fSamplingTime;      // time step in ns
  Double_t fT0;                // T0
  Double_t fDf;                // Frequency
  Double_t fA0;                // A0
  Double_t *fZeroM;            // ZeroM
  Double_t *fZeroR;            // ZeroR
  Double_t *fZeroI;            // ZeroI
  Double_t *fPoleM;            // PoleM
  Double_t *fPoleR;            // PoleR
  Double_t *fPoleI;            // PoleI
  Double_t *fTfR;     // Transfer function (real part)
  Double_t *fTfI;     // Transfer function (imaginary part)
  Double_t *fWR;     // Fourier Weights (real part)
  Double_t *fWI;     // Fourier Weights (imaginary part)
  
  ClassDef(AliITSetfSDD,1)  // Class for SDD electornics
    };
    
#endif
  

