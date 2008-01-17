#ifndef ALIMUONVTRACKERDATA_H
#define ALIMUONVTRACKERDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMUONVTrackerData
/// \brief Base class for MUON data that can be presented at different levels
/// in the hierarchy of the MUON system.
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif
#ifndef ROOT_RQ_OBJECT
#   include <RQ_OBJECT.h>
#endif

class AliMUONVStore;

class AliMUONVTrackerData : public TNamed
{
  RQ_OBJECT("AliMUONVTrackerData")
  
public:
  
  AliMUONVTrackerData(const char* name="",const char* title="", Bool_t runnable=kTRUE);
  virtual ~AliMUONVTrackerData();
  
  /// Add values for one full store
  virtual Bool_t Add(const AliMUONVStore& store) = 0;

  virtual Double_t BusPatch(Int_t busPatchId, Int_t dim=0) const = 0;
  
  virtual Double_t Chamber(Int_t chamberId, Int_t dim=0) const = 0;
  
  virtual Double_t Channel(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                           Int_t dim=0) const = 0;
  
  virtual void Clear(Option_t* opt="") = 0;
  
  virtual Double_t Count(Int_t detElemId, Int_t manuId, Int_t manuChannel) const = 0;

  virtual Double_t DetectionElement(Int_t detElemId, Int_t dim=0) const = 0;
  
  virtual TString DimensionName(Int_t dim) const = 0;

  virtual Bool_t HasBusPatch(Int_t busPatchId) const = 0;

  virtual Bool_t HasChamber(Int_t chamberId) const = 0;
  
  virtual Bool_t HasDetectionElement(Int_t detElemId) const = 0;
  
  virtual Bool_t HasManu(Int_t detElemId, Int_t manuId) const = 0;

  virtual Bool_t HasPCB(Int_t detElemId, Int_t pcbIndex) const = 0;
  
  virtual Bool_t IsRunnable() const = 0;
  
  virtual Double_t Manu(Int_t detElemId, Int_t manuId, Int_t dim=0) const = 0;
  
  virtual Int_t NumberOfDimensions() const = 0;

  virtual Int_t NumberOfEvents() const = 0;

  virtual void NumberOfEventsChanged(); // *SIGNAL*
  
  const char* Name() const { return GetName(); }
  
  virtual Double_t PCB(Int_t detElemId, Int_t pcbIndex, Int_t dim=0) const = 0;
  
  /// Print all objects whose name matches wildcard
  virtual void Print(Option_t* wildcard="") const;
  
  /// Print, with option, all objects whose name matches wildcard
  virtual void Print(Option_t* wildcard, Option_t* opt) const = 0;

  virtual void SetDimensionName(Int_t index, const char* value) = 0;
  
private:
  /// not implemented
  AliMUONVTrackerData(const AliMUONVTrackerData& rhs);
  /// not implemented
  AliMUONVTrackerData& operator=(const AliMUONVTrackerData& rhs);
  
  ClassDef(AliMUONVTrackerData,1) // Base class of MUON data that can be represented graphically
};

#endif
