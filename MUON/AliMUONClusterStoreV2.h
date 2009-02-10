#ifndef ALIMUONCLUSTERSTOREV2_H
#define ALIMUONCLUSTERSTOREV2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONClusterStoreV2
/// \brief Implementation of VClusterStore
/// 
// Author Philippe Pillot, Subatech

#ifndef ALIMUONVCLUSTERSTORE_H
#  include "AliMUONVClusterStore.h"
#endif

#include "AliMUONVCluster.h"
#include <TClonesArray.h>

class AliMUONClusterStoreV2 : public AliMUONVClusterStore
{
  friend class AliMUONClusterStoreV2Iterator;
  
public:
  AliMUONClusterStoreV2();
  AliMUONClusterStoreV2(const AliMUONClusterStoreV2& store);
  AliMUONClusterStoreV2& operator=(const AliMUONClusterStoreV2& store);  
  virtual ~AliMUONClusterStoreV2();
  
  virtual void Clear(Option_t* opt="");
  
  /// Whether the Connect(TTree&) method is implemented
  virtual Bool_t CanConnect() const { return kTRUE; }
  virtual Bool_t Connect(TTree& tree, Bool_t alone=kTRUE) const;
  
  /// Create an empty copy of this
  virtual AliMUONClusterStoreV2* Create() const { return new AliMUONClusterStoreV2; }
  
  virtual AliMUONVCluster* CreateCluster(Int_t chamberId, Int_t detElemId, Int_t clusterIndex) const;
  
  using AliMUONVClusterStore::Add;
  
  virtual AliMUONVCluster* Add(const AliMUONVCluster& Cluster);
  virtual AliMUONVCluster* Add(Int_t chamberId, Int_t detElemId, Int_t clusterIndex);

  virtual AliMUONVCluster* Remove(AliMUONVCluster& cluster);

  using AliMUONVClusterStore::GetSize;
  
  /// Return the number of clusters we hold
  virtual Int_t GetSize() const {return fClusters->GetLast()+1;}
  
  using AliMUONVStore::FindObject;
  
  AliMUONVCluster* FindObject(const TObject* object) const;
  AliMUONVCluster* FindObject(UInt_t uniqueID) const;
  
  virtual TIterator* CreateIterator() const;
  virtual TIterator* CreateChamberIterator(Int_t firstChamberId, Int_t lastChamberId) const;
  
private:
  void ReMap();
  void UpdateMap(AliMUONVCluster& cluster);
  
private:
  TClonesArray* fClusters; ///< collection of clusters
  TClonesArray* fMap;      //!< index map for fast cluster retrieval
  Bool_t        fMapped;   //!< whether our internal indices are uptodate
  
  ClassDef(AliMUONClusterStoreV2,1) // Implementation of VClusterStore
};

#endif
