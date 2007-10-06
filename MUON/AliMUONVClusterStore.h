#ifndef ALIMUONVCLUSTERSTORE_H
#define ALIMUONVCLUSTERSTORE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVClusterStore
/// \brief Interface of a cluster container
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVSTORE_H
#  include "AliMUONVStore.h"
#endif

#ifndef ALIMUONVCLUSTER_H
#  include "AliMUONVCluster.h" // must be there for covariant return type of FindObjet methods
#endif

class AliMUONVCluster;

class AliMUONVClusterStore : public AliMUONVStore
{
public:
  AliMUONVClusterStore();
  virtual ~AliMUONVClusterStore();
  
  virtual Bool_t Add(TObject* object);

  /// Add a cluster object to the store
  virtual Bool_t Add(const AliMUONVCluster& Cluster) = 0;
  /// Create a new cluster with an unique ID and add it to the store
  virtual AliMUONVCluster* Add(Int_t chamberId, Int_t detElemId, Int_t clusterIndex) = 0;

  using AliMUONVStore::Create;
  
  static AliMUONVClusterStore* Create(TTree& tree);
  
  /// Create a cluster
  virtual AliMUONVCluster* CreateCluster(Int_t chamberId, Int_t detElemId, Int_t clusterIndex) const = 0;
  
  /// Return an iterator to loop over the whole store
  virtual TIterator* CreateIterator() const = 0;

  /// Return an iterator to loop over the store in the given chamber range
  virtual TIterator* CreateChamberIterator(Int_t firstChamberId, Int_t lastChamberId) const = 0;
  
  /// Clear container
  virtual void Clear(Option_t* opt="") = 0;
  
  /// Remove a cluster object to the store
  virtual AliMUONVCluster* Remove(AliMUONVCluster& cluster) = 0;
    
  using AliMUONVStore::FindObject;

  // Find an object (default is the same as in AliMUONVStore)
  virtual AliMUONVCluster* FindObject(const TObject* object) const;
  
  // Find an object by its uniqueID (default is the same as in AliMUONVStore)
  virtual AliMUONVCluster* FindObject(UInt_t uniqueID) const;
  
  ClassDef(AliMUONVClusterStore,1) // Cluster container interface
};

#endif

