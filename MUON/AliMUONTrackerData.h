#ifndef ALIMUONTRACKERDATA_H
#define ALIMUONTRACKERDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONTrackerData
/// \brief Implementation of AliMUONVTrackerData
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVTRACKERDATA_H
#  include "AliMUONVTrackerData.h"
#endif

class AliMUONVCalibParam;
class AliMUONVStore;

class AliMUONTrackerData : public AliMUONVTrackerData
{
public:
  AliMUONTrackerData(const char* name="", const char* title="", 
                     Int_t dimension=0,
                     Bool_t runnable=kTRUE);
  virtual ~AliMUONTrackerData();

  virtual Bool_t Add(const AliMUONVStore& channelValues);
  
  virtual Double_t BusPatch(Int_t busPatchId, Int_t dim=0) const;

  virtual Double_t Chamber(Int_t chamberId, Int_t dim=0) const;

  virtual Double_t Channel(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                           Int_t dim=0) const;
  
  virtual void Clear(Option_t* opt="");
  
  virtual Double_t Count(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;
  
  virtual Double_t DetectionElement(Int_t detElemId, Int_t dim=0) const;

  virtual TString DimensionName(Int_t dim) const;

  virtual Bool_t HasChamber(Int_t chamberId) const;
  
  virtual Bool_t HasBusPatch(Int_t busPatchId) const;

  virtual Bool_t HasDetectionElement(Int_t detElemId) const;

  virtual Bool_t HasManu(Int_t detElemId, Int_t manuId) const;

  virtual Bool_t HasPCB(Int_t detElemId, Int_t pcbIndex) const;
  
  /// Whether we can be run
  virtual Bool_t IsRunnable() const { return fIsRunnable; }
  
  virtual Double_t Manu(Int_t detElemId, Int_t manuId, Int_t dim=0) const;
      
  /// Returns the number of dimensions (i.e. the number of values) each element has
  virtual Int_t NumberOfDimensions() const;
  
  /// Returns the number of events we have seen so far
  virtual Int_t NumberOfEvents() const { return fNevents; }
  
  virtual Double_t PCB(Int_t detElemId, Int_t pcbIndex, Int_t dim=0) const;

  using TObject::Print;
  
  /// Print, with option, all objects whose name matches wildcard
  virtual void Print(Option_t* wildcard, Option_t* opt) const;
  
  virtual void SetDimensionName(Int_t index, const char* value);  

//protected: FIXME: uncomment when debug done

  virtual AliMUONVCalibParam* BusPatchParam(Int_t busPatch) const;

  virtual AliMUONVCalibParam* ChamberParam(Int_t chamberId) const;
  
  virtual AliMUONVCalibParam* ChannelParam(Int_t detElemId, Int_t manuId) const;

  virtual AliMUONVCalibParam* DetectionElementParam(Int_t detElemId) const;
  
  virtual AliMUONVCalibParam* ManuParam(Int_t detElemId, Int_t manuId) const;

  virtual AliMUONVCalibParam* PCBParam(Int_t detElemId, Int_t pcbIndex) const;
  
  /// Index of the dimension containing the number of time an item was hit
  virtual Int_t IndexOfNumberDimension() const { return fDimension - 1; }

  /// Index of the dimension containing the occupancy number
  virtual Int_t IndexOfOccupancyDimension() const { return fDimension - 2; }

private:
  /// Not implemented
  AliMUONTrackerData(const AliMUONTrackerData& rhs);
  /// Not implemented
  AliMUONTrackerData& operator=(const AliMUONTrackerData& rhs);
  
  AliMUONVCalibParam* CreateDouble(const AliMUONVCalibParam& param) const;

  /// Convert from external to internal index
  Int_t External2Internal(Int_t index) const { return index*2; }

  Bool_t InternalAdd(const AliMUONVStore& channelValues);

  void SetInternalDimensionName(Int_t index, const char* value);  

  Double_t Value(const AliMUONVCalibParam& param, Int_t i, Int_t dim) const;
    
private:
    
  AliMUONVStore* fChannelValues; ///< the channel store
  AliMUONVStore* fManuValues; ///< the manu store
  AliMUONVStore* fBusPatchValues; ///< the bus patch store
  AliMUONVStore* fDEValues; ///< the detection element store
  AliMUONVStore* fChamberValues; ///< the chamber store
  AliMUONVStore* fPCBValues; ///< the pcb store
  Int_t fDimension; ///< the dimension of the data
  Int_t fNevents; ///< the number of events treated
  TObjArray* fDimensionNames; ///< the names of the dimensions
  Int_t fExternalDimension; ///< number of interface values per item 
  Bool_t fIsRunnable; ///< whether we can deal with more than one event
  
  static const Int_t fgkExtraDimension; ///< to hold extra information
  static const Int_t fgkVirtualExtraDimension; ///< to give access to information not stored, but computed on the fly
  
  ClassDef(AliMUONTrackerData,1) // Implementation of AliMUONVTrackerData
};

#endif
