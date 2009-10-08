#ifndef ALIMUONQAMAPPINGCHECK_H
#define ALIMUONQAMAPPINGCHECK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONQAMappingCheck
/// \brief Class to be called from AliMUONQADataMakerRec
/// 
// author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONDigitCalibrator;
class AliMUONGeometryTransformer;
class AliMUONVCluster;
class AliMUONVStore;
class AliMUONVTrackerData;

class AliMUONQAMappingCheck : public TObject
{
public:
  AliMUONQAMappingCheck(Int_t runNumber);
  virtual ~AliMUONQAMappingCheck();

  void NewEvent();
  
  AliMUONVTrackerData* CreateData(const char* name) const;

  void Store(AliMUONVCluster& cluster);
  
private:
  
  /// not defined on purpose
  AliMUONQAMappingCheck(const AliMUONQAMappingCheck& rhs);
  /// not defined on purpose
  AliMUONQAMappingCheck& operator=(const AliMUONQAMappingCheck& rhs);
  
  void GetClusterLocation(AliMUONVCluster& cluster, 
                          Int_t& manuBending, Int_t& manuBendingChannel, 
                          Int_t& manuNonBending, Int_t& manuNonBendingChannel,
                          Bool_t& monoCathode, Bool_t& legitimateMonoCathode);

  void AddClusterLocation(Int_t detElemId,
                          Int_t manuId, Int_t manuChannel, 
                          Bool_t monoCathode, Bool_t legitimateMonoCathode);
  
  Bool_t IsChannelDead(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;
  
  Bool_t IsManuDead(Int_t detElemId, Int_t manuId) const;
  
private:
   
  AliMUONVStore* fStore; //!< store cluster informations at manu level
  
  AliMUONGeometryTransformer* fGeometryTransformer;  //!< to go from global to local DE coordinates
  
  AliMUONDigitCalibrator* fDigitCalibrator; //!< to get statusmap
  
  Int_t fNumberOfEvents; //!< number of events seen
  
  Int_t fNumberOfClusters; //!< total number of clusters seen
  
  Int_t fNumberOfMonoCathodeClusters; //!< total number of mono-cathode clusters seen
  
  Int_t fNumberOfLegitimateMonoCathodeClusters; //!< number of mono-cathode that have a reason to be so
  
  ClassDef(AliMUONQAMappingCheck,1) // QADataMaker helper class
};

#endif
