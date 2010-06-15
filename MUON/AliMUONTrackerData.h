#ifndef ALIMUONTRACKERDATA_H
#define ALIMUONTRACKERDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONTrackerData
/// \brief Implementation of AliMUONVTrackerData
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVTRACKERDATA_H
#  include "AliMUONVTrackerData.h"
#endif

class AliMUONSparseHisto;
class AliMUONVCalibParam;
class AliMUONVStore;
class AliMpDetElement;

class AliMUONTrackerData : public AliMUONVTrackerData
{
public:
  AliMUONTrackerData(const char* name="", const char* title="", 
                     Int_t dimension=0,
                     Bool_t issingleevent=kFALSE);
  
  AliMUONTrackerData(const char* name, const char* title,
                     const AliMUONVStore& manuValues);
  
  virtual ~AliMUONTrackerData();

  Bool_t Add(const AliMUONTrackerData& data);
  
  virtual Bool_t Add(const AliMUONVStore& channelValues, TArrayI* nofEventsPerDDL=0x0);

  virtual Bool_t Replace(const AliMUONVStore& channelValues);

  virtual Double_t BusPatch(Int_t busPatchId, Int_t dim=0) const;

  virtual Double_t Chamber(Int_t chamberId, Int_t dim=0) const;

  virtual Double_t Channel(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                           Int_t dim=0) const;
  
  virtual void Clear(Option_t* opt="");
  
  virtual Double_t Count(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;
  
  virtual Double_t DetectionElement(Int_t detElemId, Int_t dim=0) const;

  virtual TString DimensionName(Int_t dim) const;
  
  virtual TString ExternalDimensionName(Int_t dim) const;

  virtual Bool_t HasChamber(Int_t chamberId) const;
  
  virtual Bool_t HasBusPatch(Int_t busPatchId) const;

  virtual Bool_t HasDetectionElement(Int_t detElemId) const;

  virtual Bool_t HasManu(Int_t detElemId, Int_t manuId) const;

  virtual Bool_t HasPCB(Int_t detElemId, Int_t pcbIndex) const;
  
  /// Whether we can be run
  virtual Bool_t IsSingleEvent() const { return fIsSingleEvent; }
  
  virtual Double_t Manu(Int_t detElemId, Int_t manuId, Int_t dim=0) const;
      
  /// Returns the number of dimensions (i.e. the number of values) each element has
  virtual Int_t NumberOfDimensions() const;
  
  /// The number of values we are inputting
  virtual Int_t ExternalDimension() const { return fExternalDimension; }

  /// Convert from internal to external dimension
  virtual Int_t InternalToExternal(Int_t dim) const { return dim/2; }

  /// Returns the number of events we have seen so far
  virtual Int_t NumberOfEvents(Int_t ddlNumber) const;
  
  virtual Double_t PCB(Int_t detElemId, Int_t pcbIndex, Int_t dim=0) const;

  using TObject::Print;
  
  /// Print, with option, all objects whose name matches wildcard
  virtual void Print(Option_t* wildcard, Option_t* opt) const;
  
  virtual void SetDimensionName(Int_t index, const char* value);  

  Bool_t CanHistogram() const { return kTRUE; }
  
  void MakeHistogramForDimension(Int_t index, Bool_t value, Double_t xmin=0.0, Double_t xmax=4096.0);
  
  virtual void HistogramRange(Double_t& xmin, Double_t& xmax) const { xmin = fXmin; xmax = fXmax; }

  AliMUONSparseHisto* GetManuSparseHisto(Int_t detElemId, Int_t manuId, 
                                         Int_t dim=0);

  AliMUONSparseHisto* GetManuSparseHisto(Int_t detElemId, Int_t manuId, 
                                         Int_t dim=0) const;
  
  AliMUONSparseHisto* GetChannelSparseHisto(Int_t detElemId, Int_t manuId, 
                                            Int_t manuChannel, Int_t dim=0);
  
  virtual AliMUONSparseHisto* GetChannelSparseHisto(Int_t detElemId, Int_t manuId, 
                                                    Int_t manuChannel, Int_t dim=0) const;

	/// Disable storing values at the channel level
	virtual void DisableChannelLevel();
	
	/// Whether we store values at the channel level or not
	virtual Bool_t IsChannelLevelEnabled() const { return fIsChannelLevelEnabled; }

  /// Disable storing values at the manu level
	virtual void DisableManuLevel();
	
	/// Whether we store values at the manu level or not
	virtual Bool_t IsManuLevelEnabled() const { return fIsManuLevelEnabled; }
  
  /// To allow merging of different objects
  virtual Long64_t Merge(TCollection* list);
    
  Bool_t ExportAsASCIIOccupancyFile(const char* filename, Int_t runNumber) const;
  
private:
    
  void FillHisto(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                 Int_t dim, Double_t value);

  AliMUONVCalibParam* BusPatchParam(Int_t busPatch, Bool_t create=kFALSE) const;

  AliMUONVCalibParam* CreateBusPatchParam(Int_t busPatch) const;
  
  AliMUONVCalibParam* ChamberParam(Int_t chamberId, Bool_t create=kFALSE) const;

  AliMUONVCalibParam* CreateChamberParam(Int_t chamberId) const;
  
  AliMUONVCalibParam* ChannelParam(Int_t detElemId, Int_t manuId,
                                   const AliMUONVCalibParam* external=0x0) const;

  AliMUONVCalibParam* DetectionElementParam(Int_t detElemId, Bool_t create=kFALSE) const;

  AliMUONVCalibParam* CreateDetectionElementParam(Int_t detElemId) const;
  
  AliMUONVCalibParam* ManuParam(Int_t detElemId, Int_t manuId, Bool_t create=kFALSE) const;

  AliMUONVCalibParam* CreateManuParam(Int_t detElemInd, Int_t manuId) const;
  
  AliMUONVCalibParam* PCBParam(Int_t detElemId, Int_t pcbIndex, Bool_t create=kFALSE) const;

  AliMUONVCalibParam* CreatePCBParam(Int_t detElemId, Int_t pcbIndex) const;
  
  /// Index of the dimension containing the number of time an item was hit
  virtual Int_t IndexOfNumberDimension() const { return fDimension - 1; }

  /// Index of the dimension containing the occupancy number
  virtual Int_t IndexOfOccupancyDimension() const { return fDimension - 2; }

  /// Whether we have histograms for a given dimension, or not
  virtual Bool_t IsHistogrammed(Int_t dim) const { return ( fHistogramming[dim] > 0 ); }

  Int_t DdlIdFromBusPatchId(Int_t buspatchid) const;
  Int_t DdlIdFromDetElemId(Int_t detelemid) const;
  Int_t DdlIdFromChamberId(Int_t chamberid) const;
  
  /// Not implemented
  AliMUONTrackerData(const AliMUONTrackerData& rhs);
  /// Not implemented
  AliMUONTrackerData& operator=(const AliMUONTrackerData& rhs);
  
  AliMUONVCalibParam* CreateDouble(const AliMUONVCalibParam& param, Int_t detElemId, Int_t manuId) const;

  Int_t GetParts(AliMUONVCalibParam* external,
                 AliMUONVCalibParam*& chamber,
                 AliMUONVCalibParam*& de,
                 AliMUONVCalibParam*& busPatch,
                 AliMUONVCalibParam*& pcb,
                 AliMUONVCalibParam*& manu,
                 AliMUONVCalibParam*& channel,
                 AliMpDetElement*& mpde);

  /// Convert from external to internal index
  Int_t External2Internal(Int_t index) const;

  void SetInternalDimensionName(Int_t index, const char* value);  

  void SetExternalDimensionName(Int_t index, const char* value);  

  Double_t Value(const AliMUONVCalibParam& param, Int_t i, Int_t dim, Int_t ddlId) const;
  
  /// The number of values we actually *store* for each item
  Int_t Dimension() const { return fDimension; }
    
  Bool_t InternalAdd(const AliMUONVStore& store, TArrayI* nevents, Bool_t replace);

  void GetDEManu(const AliMUONVCalibParam& param,
                  Int_t& detElemId, Int_t& manuId) const;
  
  void AddCalibParams(const AliMUONVCalibParam& src, AliMUONVCalibParam& dest) const;

  void Add2D(const AliMUONVStore& src, AliMUONVStore& dest) const;
  
  void Add1D(const AliMUONVStore& src, AliMUONVStore& dest) const;
  
  void AssertStores();
  
  Bool_t UpdateNumberOfEvents(TArrayI* nevents);
  
private:
  
  Bool_t fIsSingleEvent; ///< whether we can deal with more than one event
  AliMUONVStore* fChannelValues; ///< the channel store
  AliMUONVStore* fManuValues; ///< the manu store
  AliMUONVStore* fBusPatchValues; ///< the bus patch store
  AliMUONVStore* fDEValues; ///< the detection element store
  AliMUONVStore* fChamberValues; ///< the chamber store
  AliMUONVStore* fPCBValues; ///< the pcb store
  Int_t fDimension; ///< the dimension of the data
  Int_t fNevents; ///< the number of events treated
  TObjArray* fDimensionNames; ///< the names of the (internal) dimensions
  TObjArray* fExternalDimensionNames; ///< the names of the external (i.e. original) dimensions
  Int_t fExternalDimension; ///< number of interface values per item 
  /// whether we should histogram the dimension(s)
  Int_t* fHistogramming; //[fExternalDimension] whether we should histogram the dimension(s)
  AliMUONVStore* fHistos; ///< the lowest histograms we have
  Double_t fXmin; ///< min x value for histograms
  Double_t fXmax; ///< max x value for histograms
  static const Int_t fgkExtraDimension; ///< to hold extra information
  static const Int_t fgkVirtualExtraDimension; ///< to give access to information not stored, but computed on the fly
  Bool_t fIsChannelLevelEnabled; ///< whether we allow storing of channel (fChannelValues) values
  Bool_t fIsManuLevelEnabled; ///< whether we allow storing of manu (fManuValues) values
  Int_t fNofDDLs; ///< nof of DDLs we're dealing with
  /// the number of events treated (per DDL)
  Int_t* fNofEventsPerDDL; //[fNofDDLs] the number of events treated (per DDL)

  ClassDef(AliMUONTrackerData,7) // Implementation of AliMUONVTrackerData
};

#endif
