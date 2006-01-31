/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONCalibrationData
/// \brief Single entry point to access pedestals and gains from the
/// (de)calibrator or any class needing the calibration data
/// 
/// \author Laurent Aphecetche

#ifndef ALIMUONCALIBRATIONDATA_H
#define ALIMUONCALIBRATIONDATA_H

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliCDBEntry;
class AliMUONCalibParam;
class AliMUONV3DStore;

class AliMUONCalibrationData : public TObject
{
public:
  AliMUONCalibrationData(Int_t runNumber=-1, Bool_t deferredInitialization=kTRUE);
  virtual ~AliMUONCalibrationData();
  
  AliMUONCalibParam* Gain(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;
  
  Bool_t IsValid() const;
  
  AliMUONCalibParam* Pedestal(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;
  
  virtual void Print(Option_t* opt="") const;
  
  Int_t RunNumber() const;  
  
private:
  AliCDBEntry* GetEntry(const char* path) const;
  AliMUONV3DStore* Gains() const;
  AliMUONV3DStore* Pedestals() const;

private:  
  mutable Bool_t fIsValid;
  Int_t fRunNumber;
  mutable AliMUONV3DStore* fGains; //!
  mutable AliMUONV3DStore* fPedestals; //!
  
  ClassDef(AliMUONCalibrationData,1) // Storage for all MUON calibration data.
};

#endif
