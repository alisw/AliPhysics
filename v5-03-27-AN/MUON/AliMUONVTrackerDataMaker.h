#ifndef ALIMUONVTRACKERDATAMAKER_H
#define ALIMUONVTRACKERDATAMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVTrackerDataMaker
/// \brief Producer of some AliMUONVTrackerData
/// 
//  Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVTrackerData;
class TCollection;

class AliMUONVTrackerDataMaker : public TObject
{
public:
  AliMUONVTrackerDataMaker();
  virtual ~AliMUONVTrackerDataMaker();
  
  /// Whether we are valid or not
  virtual Bool_t IsValid() const = 0;
  
  /// Our data
  virtual AliMUONVTrackerData* Data() const = 0;

  /// Whether or not we're the owner of our AliMUONVTrackerData
  virtual void SetOwnerOfData(Bool_t /*flag*/) { }
  
  /// Whether we can be run
  virtual Bool_t IsRunnable() const = 0;

  /// Whether we are running (must be false if IsRunnable is false)
  virtual Bool_t IsRunning() const = 0;
  
  /// Whether we're only handling event-by-event data (i.e. no accumulation)
  virtual Bool_t IsEventByEvent() const { return kFALSE; }
  
  /// Set event-by-event mode
  virtual void SetEventByEvent(Bool_t /*flag*/) { }
  
  /// Set the running state (no effect if not runnable)
  virtual void SetRunning(Bool_t flag) = 0;
  
	/// Process current event
	virtual Bool_t ProcessEvent() = 0;
	
  /// Advance to next event and process it (no effect if not runnable)
  virtual Bool_t NextEvent() { return ProcessEvent(); }
  
  /// Rewind events (no effect if not runnable)
  virtual void Rewind() = 0;
  
  /// Set our source URI
  virtual void SetSource(const char* source) = 0;
  
  /// Get our source URI
  virtual TString Source() const = 0;

  /// Get the number of events we have seen (but not necessarily used...)
  virtual Int_t NumberOfEvents() const = 0;
  
  /// Merge
  virtual Long64_t Merge(TCollection* list) = 0;
  
  /// Set event range (if not event by event)
  virtual void SetEventRange(Int_t /* firstevent */, Int_t /* lastevent */) {}
  
  ClassDef(AliMUONVTrackerDataMaker,1) // Producer of AliMUONVTrackerData
};

#endif
