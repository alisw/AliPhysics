#ifndef ALIMUONTRACKERRAWDATAMAKER_H
#define ALIMUONTRACKERRAWDATAMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONTrackerRawDataMaker
/// \brief Creator of AliMUONVTrackerData from AliRawReader
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVTRACKERDATAMAKER_H
#  include "AliMUONVTrackerDataMaker.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliRawReader;
class AliMUONCalibrationData;
class AliMUONDigitCalibrator;
class AliMUONDigitMaker;
class AliMUONVTrackerData;
class AliMUONVStore;
class AliMUONVDigitStore;

class AliMUONTrackerRawDataMaker : public AliMUONVTrackerDataMaker
{
public:
  AliMUONTrackerRawDataMaker(AliRawReader* reader = 0x0, const char* cdbpath=0x0);
  virtual ~AliMUONTrackerRawDataMaker();
  
  /// Whether we have a valid raw reader
  Bool_t IsValid() const { return fRawReader != 0x0; }
  
  /// Our data
  AliMUONVTrackerData* Data() const { return fAccumulatedData; }
  
  /// We can be run
  virtual Bool_t IsRunnable() const { return kTRUE; }
  
  /// Whether we are running or not
  virtual Bool_t IsRunning() const { return fIsRunning; }
  
  /// Set the running status
  virtual void SetRunning(Bool_t flag) { fIsRunning = flag; }
  
  Bool_t NextEvent();
  
  void Print(Option_t* opt="") const;
  
  void Rewind();
  
  /// Tell if we are owner of our data or not
  void SetOwner(Bool_t flag) { fIsOwner = flag; }
  
  /// Get our source URI
  virtual TString Source() const { return fSource.Data(); }
  
  /// Set our source URI
  void SetSource(const char* source) { fSource = source; }
  
  /// Get our digit store
  AliMUONVDigitStore* DigitStore() const { return fDigitStore; }
  
private:
  /// Not implemented
  AliMUONTrackerRawDataMaker(const AliMUONTrackerRawDataMaker& rhs);
  /// Not implemented
  AliMUONTrackerRawDataMaker& operator=(const AliMUONTrackerRawDataMaker& rhs);
  
  void ConvertDigits();
  
private:
  AliRawReader* fRawReader; ///< reader of the data (owner)
  AliMUONVTrackerData* fAccumulatedData; ///< data (owner if fIsOwner==kTRUE)
  AliMUONVStore* fOneEventData; ///< data for one event (owner)
  Bool_t fIsOwner; ///< whether we are owner of our data or not
  TString fSource; ///< where the data comes from
  Bool_t fIsRunning; ///< whether we are running or are paused
  AliMUONDigitMaker* fDigitMaker; ///< digit maker
  AliMUONDigitCalibrator* fDigitCalibrator; ///< digit calibrator (if calibrating data)
  AliMUONCalibrationData* fCalibrationData; ///< calibration data (if calibrating data)  
  AliMUONVDigitStore* fDigitStore; ///< digit store (if calibrating data)
  TString fCDBPath; ///< OCDB path (if calibrating data)
  static Int_t fgkCounter; ///< to count the number of instances
  
  ClassDef(AliMUONTrackerRawDataMaker,1) // Producer of AliMUONVTrackerData from raw data
};

#endif
