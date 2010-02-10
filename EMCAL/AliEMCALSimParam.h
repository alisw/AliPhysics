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
	Int_t   GetDigitThreshold()     const { return fDigitThreshold;}
	Float_t GetPinNoise()           const { return fPinNoise;}
	Float_t GetTimeResolution()     const { return fTimeResolution ; }
	Int_t   GetNADCEC()             const { return fNADCEC ; }
	Int_t   GetMeanPhotonElectron() const { return fMeanPhotonElectron ; }
	void    SetDigitThreshold(Int_t val)    { fDigitThreshold     = val ; }
	void    SetPinNoise(Float_t val)        { fPinNoise           = val ; }
	void    SetTimeResolution(Float_t val)  { fTimeResolution     = val ; }
	void    SetNADCED(Int_t val)            { fNADCEC             = val ; }
 	void    SetMeanPhotonElectron(Int_t val){ fMeanPhotonElectron = val ; }

	//Parameters used in SDigitizer
	Float_t GetA()                  const { return fA ; }
	Float_t GetB()                  const { return fB ; }
	Float_t GetECPrimaryThreshold() const { return fECPrimThreshold ; }
	void    SetA(Float_t val)                  { fA               = val ; }
	void    SetB(Float_t val)                  { fB               = val ; }
	void    SetECPrimaryThreshold(Float_t val) { fECPrimThreshold = val ;}

private:

  AliEMCALSimParam(Int_t i); //True constructor which should be called by GetInstance()

private:


  static AliEMCALSimParam * fgSimParam ; // pointer to the unique instance of the class

	// Digitizer
	Int_t   fDigitThreshold  ;      // Threshold for storing digits in EMC
	Int_t   fMeanPhotonElectron ;   // number of photon electrons per GeV deposited energy 
	Float_t fPinNoise ;             // Electronics noise in EMC
	Float_t fTimeResolution ;       // Time resolution of FEE electronics
	//Float_t fTimeThreshold ;        // Threshold to start timing for given crystall
	//Float_t fTimeSignalLength ;     // Length of the timing signal 
	Int_t   fNADCEC ;               // number of channels in EC section ADC
	
	// SDigitizer
	Float_t fA ;                     // Pedestal parameter
	Float_t fB ;                     // Slope Digitizition parameters
	Float_t fECPrimThreshold ;       // To store primary if EC Shower Elos > threshold
		
  ClassDef(AliEMCALSimParam,1)
};

#endif

