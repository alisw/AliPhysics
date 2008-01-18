#ifndef ALIMUONTRACKEROCDBDATAMAKER_H
#define ALIMUONTRACKEROCDBDATAMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONTrackerOCDBDataMaker
/// \brief
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVTrackerDataMaker_H
#  include "AliMUONVTrackerDataMaker.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONTrackerOCDBDataMaker : public AliMUONVTrackerDataMaker
{
public:
  AliMUONTrackerOCDBDataMaker(const char* ocdbPath="",
                               Int_t runNumber=0,
                               const char* type="");
  virtual ~AliMUONTrackerOCDBDataMaker();
  
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
  
  /// Whether we're owner of our data
  virtual void SetOwner(Bool_t flag) { fIsOwner = flag; }
  
  /// Set our source URI
  virtual void SetSource(const char* source) { fSource = source; }
  
  /// Get our source URI
  virtual TString Source() const { return fSource; }
  
private:
  /// Not implemented
  AliMUONTrackerOCDBDataMaker(const AliMUONTrackerOCDBDataMaker& rhs);
  /// Not implemented
  AliMUONTrackerOCDBDataMaker& operator=(const AliMUONTrackerOCDBDataMaker& rhs);
  
private:
  Bool_t fIsValid; ///< whether we have valid data
  Bool_t fIsOwner; ///< whether or not we're the owner of our data
  AliMUONVTrackerData* fData; ///< our data
  TString fSource; ///< our source
  
  ClassDef(AliMUONTrackerOCDBDataMaker,1) // Producer of AliMUONVTrackerData from OCDB
};

#endif
