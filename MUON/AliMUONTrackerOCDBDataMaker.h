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
  
  Bool_t IsValid() const { return fIsValid; }
  
  virtual AliMUONVTrackerData* Data() const { return fData; }
  
  virtual Bool_t IsRunnable() const { return kFALSE; }
  
  virtual Bool_t IsRunning() const { return kFALSE; }
  
  virtual void SetRunning(Bool_t /*flag*/) {}
  
  virtual Bool_t NextEvent() { return kTRUE; }
  
  virtual void Rewind() { }
  
  /// Whether we're owner of our data
  virtual void SetOwner(Bool_t flag) { fIsOwner = flag; }
  
  virtual void SetSource(const char* source) { fSource = source; }
  
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
