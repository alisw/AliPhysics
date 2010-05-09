#ifndef ALIMUONTRACKERCONDITIONDATAMAKER_H
#define ALIMUONTRACKERCONDITIONDATAMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONTrackerConditionDataMaker
/// \brief Producer of AliMUONVTrackerData from OCDB or ASCII condition files
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONVTRACKERDATAMAKER_H
#  include "AliMUONVTrackerDataMaker.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONVStore;
class TMap;

class AliMUONTrackerConditionDataMaker : public AliMUONVTrackerDataMaker
{
public:
  AliMUONTrackerConditionDataMaker();
  AliMUONTrackerConditionDataMaker(Int_t runNumber, const char* ocdbPath, const char* type);
  AliMUONTrackerConditionDataMaker(const char* filename, const char* type);
  AliMUONTrackerConditionDataMaker(const char* data, const char* type, Bool_t);
  virtual ~AliMUONTrackerConditionDataMaker();
  
  static AliMUONVStore* CreateStore(Int_t runNumber, const char* source, const char* type, Int_t& startOfValidity);
  
  static AliMUONVTrackerData* CreateData(const char* type, AliMUONVStore& source, Int_t startOfValidity);
    
  /// Whether we've been properly initialized or not
  Bool_t IsValid() const { return (fData != 0x0); }
  
  /// Return our data
  virtual AliMUONVTrackerData* Data() const { return fData; }
  
  /// We are not runnable (i.e. # of event is fixed = 1)
  virtual Bool_t IsRunnable() const { return kFALSE; }
  
  /// We cannot be running as we are not runnable...
  virtual Bool_t IsRunning() const { return kFALSE; }
  
  /// N/A
  virtual void SetRunning(Bool_t /*flag*/) {}
  
  /// N/A
  virtual Bool_t ProcessEvent() { return kTRUE; }
  
  /// N/A
  virtual void Rewind() { }
  
  /// Set our source URI
  virtual void SetSource(const char* source) { fSource = source; }
  
  /// Get our source URI
  virtual TString Source() const { return fSource.Data(); }
  
  /// Number of events is always 1
  Int_t NumberOfEvents() const { return 1; }
  
  virtual Long64_t Merge(TCollection* li);

private:
  /// Not implemented
  AliMUONTrackerConditionDataMaker(const AliMUONTrackerConditionDataMaker& rhs);
  /// Not implemented
  AliMUONTrackerConditionDataMaker& operator=(const AliMUONTrackerConditionDataMaker& rhs);
  
  static AliMUONVStore* CreateHVStore(TMap& m);
  static AliMUONVStore* CreateStatusMapStore(Int_t runNumber);
  static AliMUONVStore* CreateStatusStore(Int_t runNumber);
  static AliMUONVStore* PatchGainStore(const AliMUONVStore& gains);
  static AliMUONVStore* ExpandConfig(const AliMUONVStore& config);
  
private:
  AliMUONVTrackerData* fData; ///< our data
  TString fSource; ///< source name
  
  ClassDef(AliMUONTrackerConditionDataMaker,1) // Producer of AliMUONVTrackerData from condition data (either OCDB or ascii files)
};

#endif

