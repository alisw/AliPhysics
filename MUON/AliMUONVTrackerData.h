#ifndef ALIMUONVTRACKERDATA_H
#define ALIMUONVTRACKERDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
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
#ifndef ROOT_TQObject
#   include <TQObject.h>
#endif

class AliMUONSparseHisto;
class AliMUONVStore;
class TCollection;
class TArrayI;

class AliMUONVTrackerData : public TNamed, public TQObject
{
public:
  
  AliMUONVTrackerData(const char* name="",const char* title="", 
                      Bool_t issingleevent=kFALSE);
  virtual ~AliMUONVTrackerData();
  
  /// Add values for one event from one full store
  virtual Bool_t Add(const AliMUONVStore& store, TArrayI* arrayOfNofEventsPerDDL=0x0) = 0;

  /// Replace values
  virtual Bool_t Replace(const AliMUONVStore& store) = 0;
  
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
  
  /// Get the name of a given (internal) dimension
  virtual TString DimensionName(Int_t dim) const = 0;

  /// Get the name of a given (external) dimension
  virtual TString ExternalDimensionName(Int_t dim) const = 0;

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
  
  /// Whether we deal with only one event at a time
  virtual Bool_t IsSingleEvent() const = 0;
  
  /// Get the value for a given manu and given dimension
  virtual Double_t Manu(Int_t detElemId, Int_t manuId, Int_t dim=0) const = 0;
  
  /// The number of dimensions we are handling
  virtual Int_t NumberOfDimensions() const = 0;

  /// Convert from internal to external dimension
  virtual Int_t InternalToExternal(Int_t dim) const = 0;
  
  /// The number of dimensions we are inputting
  virtual Int_t ExternalDimension() const = 0;

  /** The number of events we've seen so far in a given DDL (or any DDL if param<0)
   ddlNumber is 0..19
   */
  virtual Int_t NumberOfEvents(Int_t ddlNumber) const = 0;

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
  virtual void MakeHistogramForDimension(Int_t /* index */, Bool_t /* value */,
    Double_t /*xmin*/=0.0, Double_t /*xmax*/=4096.0) { }

  /// Get histogram range
  virtual void HistogramRange(Double_t& xmin, Double_t& xmax) const { xmin=xmax=0.0; }

  /// Whether we have histograms for a given dimension, or not
  virtual Bool_t IsHistogrammed(Int_t /*dim*/) const { return kFALSE; }

  /// Get sparse histogram for a given channel
  virtual AliMUONSparseHisto* GetChannelSparseHisto(Int_t detElemId, Int_t manuId, 
                                                    Int_t manuChannel, Int_t dim=0) const = 0;

  /// Get sparse histogram for a given manu (valid only if IsChannelLevelEnabled()==kFALSE and IsManuLevelEnabled()==kTRUE)
  virtual AliMUONSparseHisto* GetManuSparseHisto(Int_t detElemId, Int_t manuId, 
                                                 Int_t dim=0) const = 0;
  
  /// To allow merging of different objects
  virtual Long64_t Merge(TCollection* list) = 0;

	/// Disable recording of information at the channel level
	virtual void DisableChannelLevel() = 0;
	
	/// Whether we store values at the channel level
	virtual Bool_t IsChannelLevelEnabled() const = 0;

  /// Disable recording of information at the manu level (and below)
	virtual void DisableManuLevel() = 0;
	
	/// Whether we store values at the channel level
	virtual Bool_t IsManuLevelEnabled() const = 0;
  
  /// Whether we store values at the bus patch level or not
  virtual Bool_t IsBusPatchLevelEnabled() const = 0;
  
  /// Whether we store values at the PCB level or not
  virtual Bool_t IsPCBLevelEnabled() const = 0;
  
private:
  /// not implemented
  AliMUONVTrackerData(const AliMUONVTrackerData& rhs);
  /// not implemented
  AliMUONVTrackerData& operator=(const AliMUONVTrackerData& rhs);
  
  ClassDef(AliMUONVTrackerData,2) // Base class of MUON data that can be represented graphically
};

#endif
