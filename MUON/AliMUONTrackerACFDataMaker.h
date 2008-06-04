#ifndef ALIMUONTRACKERACFDATAMAKER_H
#define ALIMUONTRACKERACFDATAMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONTrackerACFDataMaker
/// \brief
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVTrackerDataMaker_H
#  include "AliMUONVTrackerDataMaker.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONTrackerACFDataMaker : public AliMUONVTrackerDataMaker
{
public:
  AliMUONTrackerACFDataMaker(const char* acfPath="",
                               const char* type="");
  virtual ~AliMUONTrackerACFDataMaker();
  
  /// Whether we've been properly initialized or not
  Bool_t IsValid() const { return fIsValid; }
  
  /// Return our data
  virtual AliMUONVTrackerData* Data() const { return fData; }
  
  /// We are not runnable (i.e. # of event is fixed = 1)
  virtual Bool_t IsRunnable() const { return kFALSE; }
  
  /// We cannot be running as we are not runnable...
  virtual Bool_t IsRunning() const { return kFALSE; }
  
  /// N/A
  virtual void SetRunning(Bool_t /*flag*/) {}
  
  /// N/A
  virtual Bool_t NextEvent() { return kTRUE; }
  
  /// N/A
  virtual void Rewind() { }
  
  /// Set our source URI
  virtual void SetSource(const char* source) { fSource = source; }
  
  /// Get our source URI
  virtual TString Source() const { return fSource; }
  
  /// Number of events is always 1
    Int_t NumberOfEvents() const { return 1; }

  virtual Long64_t Merge(TCollection* li);
  
private:
  /// Not implemented
  AliMUONTrackerACFDataMaker(const AliMUONTrackerACFDataMaker& rhs);
  /// Not implemented
  AliMUONTrackerACFDataMaker& operator=(const AliMUONTrackerACFDataMaker& rhs);
  
private:
  Bool_t fIsValid; ///< whether we have valid data
  AliMUONVTrackerData* fData; ///< our data
  TString fSource; ///< our source
  
  ClassDef(AliMUONTrackerACFDataMaker,2) // Producer of AliMUONVTrackerData from ACF
};

#endif
