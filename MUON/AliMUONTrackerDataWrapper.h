#ifndef ALIMUONTRACKERDATAWRAPPER_H
#define ALIMUONTRACKERDATAWRAPPER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONTrackerDataWrapper
/// \brief Simple wrapper of AliMUONVTrackerData (for backward compatibility)
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVTRACKERDATAMAKER_H
#  include "AliMUONVTrackerDataMaker.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONTrackerDataWrapper : public AliMUONVTrackerDataMaker
{
public:
  AliMUONTrackerDataWrapper(AliMUONVTrackerData* data=0x0);
  virtual ~AliMUONTrackerDataWrapper();
  
  /// Whether we are valid or not
  virtual Bool_t IsValid() const { return kTRUE; }
  
  /// Our data
  virtual AliMUONVTrackerData* Data() const { return fData; }
  
  /// Whether we can be run
  virtual Bool_t IsRunnable() const { return kFALSE; }
  
  /// Whether we are running (must be false if IsRunnable is false)
  virtual Bool_t IsRunning() const { return kFALSE; }
  
  /// Set the running state (no effect if not runnable)
  virtual void SetRunning(Bool_t /*flag*/) {}
  
  /// Advance to next event (no effect if not runnable)
  virtual Bool_t ProcessEvent() { return kFALSE; }
  
  /// Rewind events (no effect if not runnable)
  virtual void Rewind() { }
  
  /// Set our source URI
  virtual void SetSource(const char* /*source*/) {}
  
  /// Get our source URI
  virtual TString Source() const { return ""; }
  
  /// Get the number of events we have seen (but not necessarily used...)
  virtual Int_t NumberOfEvents() const;
  
  virtual Long64_t Merge(TCollection* li);
  
  virtual void UpdateData(AliMUONVTrackerData* data=0x0) { fData = data; }
  
private:
    /// not implemented.
    AliMUONTrackerDataWrapper(const AliMUONTrackerDataWrapper& rhs);
  /// not implemented.
  AliMUONTrackerDataWrapper& operator=(const AliMUONTrackerDataWrapper& rhs);
  
private:
    AliMUONVTrackerData* fData; ///< our data (owner)
  
  ClassDef(AliMUONTrackerDataWrapper,1) // Wrapper of AliMUONVTrackerData
};

#endif
