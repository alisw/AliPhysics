#ifndef ALIMUONVTRACKERDATA_H
#define ALIMUONVTRACKERDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
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
class TH1;

class AliMUONVTrackerData : public TNamed
{
  RQ_OBJECT("AliMUONVTrackerData")
  
public:
  
  AliMUONVTrackerData(const char* name="",const char* title="", Bool_t runnable=kTRUE);
  virtual ~AliMUONVTrackerData();
  
  /// Add values for one full store
  virtual Bool_t Add(const AliMUONVStore& store) = 0;

  /// Get the value for a given buspatch and given dimension
  virtual Double_t BusPatch(Int_t busPatchId, Int_t dim=0) const = 0;
  
  /// Get the value for a given chamber and given dimension
  virtual Double_t Chamber(Int_t chamberId, Int_t dim=0) const = 0;
  
  /// Get the value for a given channel and given dimension
  virtual Double_t Channel(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                           Int_t dim=0) const = 0;
  
  /// Reset the data
  virtual void Clear(Option_t* opt="") = 0;
  
  /// Get the number of times a given channel was hit
  virtual Double_t Count(Int_t detElemId, Int_t manuId, Int_t manuChannel) const = 0;

  /// Get the value for a given DE and given dimension
  virtual Double_t DetectionElement(Int_t detElemId, Int_t dim=0) const = 0;
  
  /// Get the name of a given dimension
  virtual TString DimensionName(Int_t dim) const = 0;

  /// Whether we have data for a given buspath
  virtual Bool_t HasBusPatch(Int_t busPatchId) const = 0;

  /// Whether we have a given channel or not
  virtual Bool_t HasChannel(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;
  
  /// Whether we have data for a given chamber
  virtual Bool_t HasChamber(Int_t chamberId) const = 0;
  
  /// Whether we have data for a given detection element
  virtual Bool_t HasDetectionElement(Int_t detElemId) const = 0;
  
  /// Whether we have data for a given manu
  virtual Bool_t HasManu(Int_t detElemId, Int_t manuId) const = 0;

  /// Whether we have data for a given PCB
  virtual Bool_t HasPCB(Int_t detElemId, Int_t pcbIndex) const = 0;
  
  /// Whether we are runnable (e.g. can handle several events)
  virtual Bool_t IsRunnable() const = 0;
  
  /// Get the value for a given manu and given dimension
  virtual Double_t Manu(Int_t detElemId, Int_t manuId, Int_t dim=0) const = 0;
  
  /// The number of dimensions we are handling
  virtual Int_t NumberOfDimensions() const = 0;

  /// The number of events we've seen so far
  virtual Int_t NumberOfEvents() const = 0;

  /// Signal to indicate that the number of events changed
  virtual void NumberOfEventsChanged(); // *SIGNAL*
  
  /// Get our name
  const char* Name() const { return GetName(); }
  
  /// Get the value for a given PCDB and given dimension
  virtual Double_t PCB(Int_t detElemId, Int_t pcbIndex, Int_t dim=0) const = 0;
  
  /// Print all objects whose name matches wildcard
  virtual void Print(Option_t* wildcard="") const;
  
  /// Print, with option, all objects whose name matches wildcard
  virtual void Print(Option_t* wildcard, Option_t* opt) const = 0;

  /// Set the name of a given dimension
  virtual void SetDimensionName(Int_t index, const char* value) = 0;

  /// Whether or not we can make histograms.
  virtual Bool_t CanHistogram() const { return kFALSE; }
  
  /// Select a dimension to be histogrammed (if CanHistogram==kTRUE) only
  virtual void SetHistogramDimension(Int_t /* index */, Bool_t /* value */) { }

  /// Create (if CanHistogram) an histo for a given channel
  virtual TH1* CreateChannelHisto(Int_t /*detElemId*/, Int_t /*manuId*/, 
                                  Int_t /*manuChannel*/, Int_t /*dim*/=0) { return 0x0; }
  
  /// Create (if CanHistogram) an histo for a given bus patch
  virtual TH1* CreateBusPatchHisto(Int_t /*busPatchId*/, Int_t /*dim*/=0)  { return 0x0; }
  
  /// Create (if CanHistogram) an histo for a given detection element
  virtual TH1* CreateDEHisto(Int_t /*detElemId*/, Int_t /*dim*/=0)  { return 0x0; }
  
  /// Create (if CanHistogram) an histo for a given manu
  virtual TH1* CreateManuHisto(Int_t /*detElemId*/, Int_t /*manuId*/, Int_t /*dim*/=0)  { return 0x0; }
  
  /// Create (if CanHistogram) an histo for a given pcb
  virtual TH1* CreatePCBHisto(Int_t /*detElemId*/, Int_t /*pcbIndex*/, Int_t /*dim*/=0)  { return 0x0; }
  
  /// Create (if CanHistogram) an histo for a given chamber
  virtual TH1* CreateChamberHisto(Int_t /*chamberId*/, Int_t /*dim*/=0)  { return 0x0; }
  
private:
  /// not implemented
  AliMUONVTrackerData(const AliMUONVTrackerData& rhs);
  /// not implemented
  AliMUONVTrackerData& operator=(const AliMUONVTrackerData& rhs);
  
  ClassDef(AliMUONVTrackerData,1) // Base class of MUON data that can be represented graphically
};

#endif
