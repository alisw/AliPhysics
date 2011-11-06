#ifndef ALIEMCALSIMPARAM_H
#define ALIEMCALSIMPARAM_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id: AliEMCALSimParam.h  $ */
/*
//
// Base class for the EMCAL simulation parameters.
//
//
*/

#include "TNamed.h"

class AliEMCALSimParam : public TNamed {

public:

  AliEMCALSimParam();
  AliEMCALSimParam(const AliEMCALSimParam& recoParam);
  AliEMCALSimParam& operator = (const AliEMCALSimParam& recoParam);
  virtual ~AliEMCALSimParam() {}

  static AliEMCALSimParam * GetInstance() ;
  virtual void Print(Option_t * option="") const ;

  //Parameters used in Digitizer
  Int_t    GetDigitThreshold()          const { return fDigitThreshold     ; }
  Float_t  GetPinNoise()                const { return fPinNoise           ; }
  Double_t GetTimeNoise()               const { return fTimeNoise          ; }
  Double_t GetTimeDelay()               const { return fTimeDelay          ; }
  Double_t GetTimeResolutionPar0()      const { return fTimeResolutionPar0 ; }
  Double_t GetTimeResolutionPar1()      const { return fTimeResolutionPar1 ; }
  Int_t    GetNADCEC()                  const { return fNADCEC             ; }
  Int_t    GetMeanPhotonElectron()      const { return fMeanPhotonElectron ; }
  Float_t  GetGainFluctuations()        const { return fGainFluctuations   ; }
  
  void     SetDigitThreshold(Int_t val)       { fDigitThreshold      = val ; }
  void     SetPinNoise(Float_t val)           { fPinNoise            = val ; }
  void     SetTimeNoise(Float_t val)          { fTimeNoise           = val ; }
  void     SetTimeDelay(Double_t val)         { fTimeDelay           = val ; }
  void     SetTimeResolutionPar0(Double_t val){ fTimeResolutionPar0  = val ; }
  void     SetTimeResolutionPar1(Double_t val){ fTimeResolutionPar1  = val ; }
  void     SetNADCED(Int_t val)               { fNADCEC              = val ; }
  void     SetMeanPhotonElectron(Int_t val)   { fMeanPhotonElectron  = val ; }
  void     SetGainFluctuations(Float_t val)   { fGainFluctuations    = val ; }

  //Parameters used in SDigitizer
  Float_t  GetA()                       const { return fA                  ; }
  Float_t  GetB()                       const { return fB                  ; }
  Float_t  GetECPrimaryThreshold()      const { return fECPrimThreshold    ; }
  
  void     SetA(Float_t val)                  { fA                   = val ; }
  void     SetB(Float_t val)                  { fB                   = val ; }
  void     SetECPrimaryThreshold(Float_t val) { fECPrimThreshold     = val ; }


private:


  static AliEMCALSimParam * fgSimParam ; // pointer to the unique instance of the class

  // Digitizer
  Int_t    fDigitThreshold  ;     // Threshold for storing digits in EMC
  Int_t    fMeanPhotonElectron ;  // number of photon electrons per GeV deposited energy 
  Float_t  fGainFluctuations ;    // correct fMeanPhotonElectron by the gain fluctuations
  Float_t  fPinNoise ;            // Electronics noise in EMC, APD
  Double_t fTimeNoise ;           // Electronics noise in EMC, time
  Double_t fTimeDelay;            // Time delay to reproduce data delay
  Double_t fTimeResolutionPar0 ;  // Time resolution of FEE electronics
  Double_t fTimeResolutionPar1 ;  // Time resolution of FEE electronics
  Int_t    fNADCEC ;              // number of channels in EC section ADC
	
  // SDigitizer
  Float_t fA ;                     // Pedestal parameter
  Float_t fB ;                     // Slope Digitizition parameters
  Float_t fECPrimThreshold ;       // To store primary if EC Shower Elos > threshold
		
  ClassDef(AliEMCALSimParam,6)
};

#endif

