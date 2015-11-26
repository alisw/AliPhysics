#ifndef ALIZDCDIGITIZER_H
#define ALIZDCDIGITIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  		Digitizer class for ZDC	      //
////////////////////////////////////////////////

#include "AliDigitizer.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliZDCPedestals.h"

class TFile;
class AliDigitizationInput;

class AliZDC;
class AliZDCHit;
class AliZDCDigit;

class AliZDCDigitizer: public AliDigitizer {

public:
  AliZDCDigitizer();
  AliZDCDigitizer(AliDigitizationInput* digInput);
  virtual ~AliZDCDigitizer();
   
  virtual Bool_t Init();
  virtual void Digitize(Option_t* option=0);    

  //  PM gain
  void    SetPMGain(Int_t det, Int_t pmDet, Int_t pmGain)
    {fPMGain[det][pmDet] = pmGain;}
  Float_t GetPMGain(Int_t det, Int_t pmDet) const
    {return fPMGain[det][pmDet];}
  //  Conversion factor from charge to ADC channels
  //	      F = 1.6E-19 / Resolution [Coulomb/ch]
  void    SetADCRes(Int_t *adcRes) {for (Int_t i=0;i<2;i++) fADCRes[i] = adcRes[i];}
  //  Two conversion factor are needed for ADC CAEN V965 
  Float_t GetADCRes(Int_t i) const {return fADCRes[i];}
  
  void	  SetCalibrationOn() {fIsCalibration=1;}  
  AliCDBStorage    *SetStorage(const char* uri);
  AliZDCPedestals  *GetPedData() const; 
  
  void    SetSpectators2Track() {fSpectators2Track=kTRUE;}
  void    SetBeamEnergy(Float_t beamEnergy) {fBeamEnergy = beamEnergy;}

  // Added for p-A simulations
  void    SetpAsystem() {fIspASystem=kTRUE;}

  // Added for RELDIS
  void    SetRELDISGenerator() {fIsRELDISgen=kTRUE;}
  
  void    SetSpectatorParam(int imod=1) {fSpectatorParam=imod;}

  void    SpectatorSignal(Int_t specType, Int_t numEvents, Float_t pm[5][5]);
  void    SpectatorsFromHijing(Int_t specType, Int_t numEvents, Float_t pm[5][5]);


private:

  AliZDCDigitizer(const AliZDCDigitizer&);
  AliZDCDigitizer& operator=(const AliZDCDigitizer&);

  void    CalculatePMTGains();
  void    ReadPMTGains();

  void    Fragmentation(Float_t impPar, Int_t specN, Int_t specP,
                        Int_t &freeSpecN, Int_t &freeSpecP) const;

  Int_t   Phe2ADCch(Int_t Detector, Int_t Quadrant, Float_t Light, 
                    Int_t Res) const;
  Int_t   Pedestal(Int_t Detector, Int_t Quadrant, Int_t Res) const;

  Float_t fPMGain[5][5];      	// PM gain
  Float_t fADCRes[2];	      	// ADC conversion factors
  Int_t   fIsCalibration; 	// !=0 if simulation creates calibration data
  Bool_t  fIsSignalInADCGate;   // true if signal in ADC gate
  Float_t fFracLostSignal;      // fraction of lost signal
  
  AliZDCPedestals  *fPedData; 	//! pedestal calibration data
  
  Bool_t  fSpectators2Track;    // should digitizer track spectators
  Float_t fBeamEnergy;          // beam energy
  TString fBeamType;		// beam type
  
  // Added for p-A simulations
  Bool_t  fIspASystem;       	// true if collision system is p-A

  // Added for RELDIS
  Bool_t  fIsRELDISgen;  	// true if generator is RELDIS
  
  // Fragmentation is derived from RUN1 data
  TFile  *fSpectatorData;	// pointer to stored spectator data files
  Int_t   fSpectatorParam;      // kinematic model fro spectators (1=from AliGenZDC =DEFAULT, 2=from HIJING)
  
  ClassDef(AliZDCDigitizer, 15)     // digitizer for ZDC
};    
#endif
