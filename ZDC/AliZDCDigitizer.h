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

  // Added for p-A simulations
  void    SetpAsystem() {fIspASystem=kTRUE;}

private:

  AliZDCDigitizer(const AliZDCDigitizer&);
  AliZDCDigitizer& operator=(const AliZDCDigitizer&);

  void    Fragmentation(Float_t impPar, Int_t specN, Int_t specP,
                        Int_t &freeSpecN, Int_t &freeSpecP) const;
  void    SpectatorSignal(Int_t SpecType, Int_t numEvents, 
                          Float_t pm[3][5]) const;

  Int_t   Phe2ADCch(Int_t Detector, Int_t Quadrant, Float_t Light, 
                    Int_t Res) const;
  Int_t   Pedestal(Int_t Detector, Int_t Quadrant, Int_t Res) const;

  Float_t fPMGain[5][5];      	// PM gain
  Float_t fADCRes[2];	      	// ADC conversion factors
  Int_t   fIsCalibration; 	// !=0 if simulation creates calibration data
  Bool_t  fIsSignalInADCGate;   // true if signal in ADC gate
  Float_t fFracLostSignal;      // fraction of lost signal
  
  AliZDCPedestals  *fPedData; 	   //! pedestal calibration data
  
  Bool_t  fSpectators2Track;    // should digitizer track spectators
  Float_t fBeamEnergy;          // beam energy taken from GRP object
  
  // Added for p-A simulations
  Bool_t  fIspASystem;       	// true if collision system is p-A
       
  ClassDef(AliZDCDigitizer, 13)     // digitizer for ZDC
};    
#endif
