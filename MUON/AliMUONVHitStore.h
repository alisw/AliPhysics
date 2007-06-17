#ifndef ALIMUONVHITSTORE_H
#define ALIMUONVHITSTORE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup sim
/// \class AliMUONVHitStore
/// \brief Virtual store to hold digit
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVSTORE_H
#  include "AliMUONVStore.h"
#endif

class AliMUONHit;
class TClonesArray;
class TCollection;

class AliMUONVHitStore : public AliMUONVStore
{
public:
  AliMUONVHitStore();
  virtual ~AliMUONVHitStore();

  /// Add an object, if of type AliMUONHit
  virtual Bool_t Add(TObject* object);
  
  /// Add a digit
  virtual void Add(const AliMUONHit& hit) = 0;
  
  using AliMUONVStore::Create;
  
  /// Create a store from the tree (if possible).
  static AliMUONVHitStore* Create(TTree& tree);
  
  /// Return an iterator to loop over hits
  virtual TIterator* CreateIterator() const = 0;
  
  /// Must be implemented to allow connection using MCApp()->AddHitList()
  virtual TCollection* Collection() = 0;
  
  using AliMUONVStore::GetSize;
  
  ClassDef(AliMUONVHitStore,1) // Base class of a MUON hit store
};

#endif
