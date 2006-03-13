/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONCalibrationData
/// \brief Single entry point to access MUON calibration data.
///
/// For the moment, this class stores pedestals, gains and deadchannels
/// that are fetched from the CDB.
/// 
/// \author Laurent Aphecetche

#ifndef ALIMUONCALIBRATIONDATA_H
#define ALIMUONCALIBRATIONDATA_H

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliCDBEntry;
class AliMUONVCalibParam;
class AliMUONV2DStore;

class AliMUONCalibrationData : public TObject
{
public:
  /** Constructor.
    * @param runNumber is used as a key to the CDB
    * @param deferredInitialization if kFALSE, all the calibrations are fetched
    * regardless of whether you'll use them or not.
    */
  AliMUONCalibrationData(Int_t runNumber=-1, Bool_t deferredInitialization=kTRUE);
  virtual ~AliMUONCalibrationData();

  /// Get the DeadChannel calibration object for channels within (detElemId,manuId).
  AliMUONVCalibParam* DeadChannel(Int_t detElemId, Int_t manuId) const;
  
  /// Get the Gain calibration object for channels within (detElemId,manuId).
  AliMUONVCalibParam* Gain(Int_t detElemId, Int_t manuId) const;

  /// Whether this object is valid or not (might be invalid if fetching from CDB failed).
  Bool_t IsValid() const { return fIsValid; }
  
  /// Get the Pedestal calibration object for channels within (detElemId,manuId).
  AliMUONVCalibParam* Pedestal(Int_t detElemId, Int_t manuId) const;
  
  /// Dump to screen.
  virtual void Print(Option_t* opt="") const;
  
  /// The runnumber used by this object.
  Int_t RunNumber() const { return fRunNumber; }
  
protected:
  AliMUONCalibrationData(const AliMUONCalibrationData& right);
  AliMUONCalibrationData&  operator = (const AliMUONCalibrationData& right);
     
private:
  AliCDBEntry* GetEntry(const char* path) const;
  AliMUONV2DStore* Gains() const;
  AliMUONV2DStore* Pedestals() const;
  AliMUONV2DStore* DeadChannels() const;
  
private:  
  mutable Bool_t fIsValid;
  Int_t fRunNumber;
  mutable AliMUONV2DStore* fGains; //!
  mutable AliMUONV2DStore* fPedestals; //!
  mutable AliMUONV2DStore* fDeadChannels; //!
  
  ClassDef(AliMUONCalibrationData,2) // Storage for all MUON calibration data.
};

#endif
