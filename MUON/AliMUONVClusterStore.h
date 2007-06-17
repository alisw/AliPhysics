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

class AliMUONRawCluster;

class AliMUONVClusterStore : public AliMUONVStore
{
public:
  AliMUONVClusterStore();
  virtual ~AliMUONVClusterStore();
  
  virtual Bool_t Add(TObject* object);

  /// Add a cluster object to the store
  virtual Bool_t Add(const AliMUONRawCluster& Cluster) = 0;

  using AliMUONVStore::Create;
  
  static AliMUONVClusterStore* Create(TTree& tree);
  
  /// Return an iterator to loop over the whole store
  virtual TIterator* CreateIterator() const = 0;

  /// Return an iterator to loop over the store in the given chamber range
  virtual TIterator* CreateChamberIterator(Int_t firstChamberId, Int_t lastChamberId) const = 0;

  /// Remove a cluster object to the store
  virtual AliMUONRawCluster* Remove(AliMUONRawCluster& cluster) = 0;
    
  using AliMUONVStore::GetSize;
  
  ClassDef(AliMUONVClusterStore,1) // Cluster container interface
};

#endif
