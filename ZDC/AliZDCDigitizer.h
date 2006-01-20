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
#include "AliZDCCalibData.h"

class AliRunDigitizer;

class AliZDC;
class AliZDCHit;
class AliZDCMergedHit;
class AliZDCDigit;

class AliZDCDigitizer: public AliDigitizer {

public:
  AliZDCDigitizer();
  AliZDCDigitizer(AliRunDigitizer* manager);
  virtual ~AliZDCDigitizer();
   
  virtual Bool_t Init();
  virtual void Exec(Option_t* option=0);    

  //  PM gain
  void    SetPMGain(Int_t det, Int_t pmDet, Int_t pmGain)
    {fPMGain[det][pmDet] = pmGain;}
  Float_t GetPMGain(Int_t det, Int_t pmDet) const
    {return fPMGain[det][pmDet];}
  //  Conversion factor from charge to ADC channels
  //	      F = 1.6E-19 / Resolution [Coulomb/ch]
  void    SetADCRes(Int_t *adcRes)
  //  Two conversion factor are needed for ADC CAEN V965 
    {for (Int_t i=0;i<2;i++) fADCRes[i] = adcRes[i];}
  Float_t GetADCRes(Int_t i) const {return fADCRes[i];}
  
  void   GetStorage(const char* uri) {fStorage = AliCDBManager::Instance()->GetStorage(uri);}
  AliCDBStorage   *SetStorage(const char* uri);
  AliZDCCalibData *GetCalibData(int runNumber) const; 

private:
  void    Fragmentation(Float_t impPar, Int_t specN, Int_t specP,
                        Int_t &freeSpecN, Int_t &freeSpecP) const;
  void    SpectatorSignal(Int_t SpecType, Int_t numEvents, 
                          Float_t pm[3][5]) const;

  Int_t   Phe2ADCch(Int_t Detector, Int_t Quadrant, Float_t Light, 
                    Int_t Res) const;
  Int_t   Pedestal(Int_t Detector, Int_t Quadrant, Int_t Res) const;

  Float_t fPMGain[3][5];      	// PM gain
  Float_t fADCRes[2];	      	// ADC conversion factors
  
  AliCDBStorage *fStorage; 	//! storage
  AliZDCCalibData *fCalibData; 	//! calibration data
       
  ClassDef(AliZDCDigitizer, 3)     // digitizer for ZDC
};    
#endif
