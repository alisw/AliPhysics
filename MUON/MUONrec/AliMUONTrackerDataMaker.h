#ifndef ALIMUONTRACKERDATAMAKER_H
#define ALIMUONTRACKERDATAMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec 
/// \class AliMUONTrackerDataMaker
/// \brief Implementation of VTrackerDataMaker to read raw data
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVTRACKERDATAMAKER_H
#  include "AliMUONVTrackerDataMaker.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONRecoParam;
class AliMUONCalibrationData;
class AliMUONDigitCalibrator;
class AliMUONVStore;
class AliMUONVTrackerData;
class AliRawReader;
class AliMUONLogger;

class AliMUONTrackerDataMaker : public AliMUONVTrackerDataMaker
{
public:
  AliMUONTrackerDataMaker(TRootIOCtor*);
  
  AliMUONTrackerDataMaker(const AliMUONRecoParam* recoParam,
                          Int_t runNumber,
                          AliRawReader* rawReader,
                          const char* cdbPath,
                          const char* calibMode,
                          Bool_t histogram=kFALSE,
                          Double_t xmin=0.0,
                          Double_t xmax=4095.0);
  
  AliMUONTrackerDataMaker(const AliMUONRecoParam* recoParam,
                          AliRawReader* rawReader,
                          const char* cdbPath,
                          const char* calibMode,
                          Bool_t histogram=kFALSE,
                          Double_t xmin=0.0,
                          Double_t xmax=4095.0);
  
  AliMUONTrackerDataMaker(AliRawReader* rawReader, Bool_t histogram=kFALSE);
  
  virtual ~AliMUONTrackerDataMaker();
  
  Bool_t Add(const AliMUONTrackerDataMaker& other);

  /// Whether we have a valid reader or not
  Bool_t IsValid() const { return fRawReader != 0x0; }
  
  /// Get our accumulated data
  AliMUONVTrackerData* Data() const { return fAccumulatedData; }
  
  /// Whether or not we're the owner of our fAccumulatedData
  void SetOwnerOfData(Bool_t flag) { fIsOwnerOfAccumulatedData = flag; }

  /// Whether we're only handling event-by-event data (i.e. no accumulation)
  Bool_t IsEventByEvent() const { return fIsEventByEvent; }
  
  /// Set event-by-event mode
  void SetEventByEvent(Bool_t flag) { fIsEventByEvent = flag; }
  
  /// We can run if we have a reader
  Bool_t IsRunnable() const { return IsValid(); }
  
  /// Whether we are running or not
  Bool_t IsRunning() const { return fIsRunning; }
  
  /// Set the runnning status
  void SetRunning(Bool_t flag) { fIsRunning = flag; }
  
  Bool_t ProcessEvent();
  
  Bool_t NextEvent();
  
  void Print(Option_t* opt="") const;
  
  void Rewind();
  
  /// Get our source URI
  TString Source() const { return fSource.Data(); }
  
  /// Set our source URI
  void SetSource(const char* source) { fSource = source; }
  
  /// Number of events seen
  Int_t NumberOfEvents() const { return fNumberOfEvents; }

  /// Number of physics events seen
  Int_t NumberOfPhysicsEvents() const { return fNumberOfPhysicsEvents; }

  /// Number of good physics events seen
  Int_t NumberOfGoodPhysicsEvents() const { return fNumberOfGoodPhysicsEvents; }

  Long64_t Merge(TCollection* li);
  
  void SetRawReader(AliRawReader* rawReader);
  
  /// Set the error logger
  void EnableErrorLogger(AliMUONLogger* logger) { fLogger = logger; }
  
  /// Whether last decoded event was empty
  Bool_t LastEventWasEmpty() const { return fLastEventWasEmpty; }
  
  /// Whether or not we should try to recover corrupted raw data
  void SetTryRecover(Bool_t flag) { fTryRecover = flag; }
  
  /// Set the event range to consider
  void SetEventRange(Int_t first, Int_t last) { fFirstEvent=first; fLastEvent=last; }

private:
  /// not implemented
  AliMUONTrackerDataMaker(const AliMUONTrackerDataMaker& rhs);
  /// not implemented
  AliMUONTrackerDataMaker& operator=(const AliMUONTrackerDataMaker& rhs);

  void Ctor(const AliMUONRecoParam* param,
            Int_t runNumber,
            const char* calibMode,
            Bool_t histogram,
            Double_t xmin=0.0, 
            Double_t xmax=4095.0);
  
private:
  AliRawReader* fRawReader; //!< reader of the data (owner or not)
  AliMUONVTrackerData* fAccumulatedData; ///< data (owner or not)
  Bool_t fIsOwnerOfAccumulatedData; ///< owner or not of fAccumulatedData
  AliMUONVStore* fOneEventData; ///< data for a single event (owner)
  AliMUONDigitCalibrator* fDigitCalibrator; //!< digit calibrator (if calibrating)
  AliMUONCalibrationData* fCalibrationData; ///< calibration data (if calibrating)
  TString fSource; ///< where the data comes from
  TString fOCDBPath; ///< OCDB path (if calibrating)
  Int_t fNumberOfEvents; ///< number of events seen
  Int_t fRunNumber; ///< run number of the data
  Bool_t fIsRunning; ///< whether we are running or not
  Bool_t fIsOwnerOfRawReader; ///< whether we must delete rawReader or not
  Bool_t fIsEventByEvent; ///< we only keep one event's data (no accumulation)
  static Int_t fgkCounter; ///< to count the number of instances
  AliMUONLogger* fLogger; ///< error logger (not owner)
  Bool_t fLastEventWasEmpty; ///< whether last decoded event was empty
  Int_t fNumberOfPhysicsEvents; ///< number of physics events seen
  Int_t fNumberOfGoodPhysicsEvents; ///< number of errors with no (fatal) readout error
  Bool_t fTryRecover; ///< whether we should try to recover corrupted raw data
  Int_t fFirstEvent; ///< first event to consider
  Int_t fLastEvent; ///< last event to consider

  ClassDef(AliMUONTrackerDataMaker,5) // Producer of AliMUONVTrackerData from raw
};

#endif
