#ifndef ALIMUONTRACKERRAWDATAMAKER_H
#define ALIMUONTRACKERRAWDATAMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONTrackerRawDataMaker
/// \brief Creator of raw AliMUONVTrackerData from AliRawReader
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVTRACKERDATAMAKER_H
#  include "AliMUONVTrackerDataMaker.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliRawReader;
class AliMUONVStore;
class AliMUONVTrackerData;

class AliMUONTrackerRawDataMaker : public AliMUONVTrackerDataMaker
{
public:
  AliMUONTrackerRawDataMaker(AliRawReader* reader = 0x0,
                             Bool_t histogram=kFALSE,
                             Bool_t useHPdecoder=kFALSE);
  virtual ~AliMUONTrackerRawDataMaker();
  
  /// Whether we have a valid raw reader
  Bool_t IsValid() const { return fRawReader != 0x0; }
  
  /// Our data
  AliMUONVTrackerData* Data() const { return fAccumulatedData; }
  
  /// Whether we're only handling event-by-event data (i.e. no accumulation)
  virtual Bool_t IsEventByEvent() const { return fIsEventByEvent; }
  
  /// Set event-by-event mode
  virtual void SetEventByEvent(Bool_t flag) { fIsEventByEvent = flag; }
  
  /// We can be run if we have a reader
  virtual Bool_t IsRunnable() const { return IsValid(); }
  
  /// Whether we are running or not
  virtual Bool_t IsRunning() const { return fIsRunning; }
  
  /// Set the running status
  virtual void SetRunning(Bool_t flag) { fIsRunning = flag; }
  
  Bool_t NextEvent();
  
  void Print(Option_t* opt="") const;
  
  void Rewind();
  
  /// Get our source URI
  virtual TString Source() const { return fSource.Data(); }
  
  /// Set our source URI
  void SetSource(const char* source) { fSource = source; }
  
  /// Number of events seen
  Int_t NumberOfEvents() const { return fNumberOfEvents; }

  Long64_t Merge(TCollection* li);
  
private:
  /// Not implemented
  AliMUONTrackerRawDataMaker(const AliMUONTrackerRawDataMaker& rhs);
  /// Not implemented
  AliMUONTrackerRawDataMaker& operator=(const AliMUONTrackerRawDataMaker& rhs);
  
private:
  AliRawReader* fRawReader; //!< reader of the data (owner)
  AliMUONVTrackerData* fAccumulatedData; ///< data (owner)
  AliMUONVStore* fOneEventData; ///< data for one event (owner)
  TString fSource; ///< where the data comes from
  Bool_t fIsRunning; ///< whether we are running or are paused
  Int_t fNumberOfEvents; ///< number of events seen
  Int_t fRunNumber; ///< run number of the data
  Bool_t fIsEventByEvent; ///< we only keep one event's data (no accumulation)
  Bool_t fUseHPDecoder; ///< whether to use high performance decoder or not (false by default)
  static Int_t fgkCounter; ///< to count the number of instances
  
  ClassDef(AliMUONTrackerRawDataMaker,3) // Producer of AliMUONVTrackerData from raw data
};

#endif
