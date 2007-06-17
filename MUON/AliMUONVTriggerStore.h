#ifndef ALIMUONVTRIGGERSTORE_H
#define ALIMUONVTRIGGERSTORE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONVTriggerStore
/// \brief Base class of a trigger information store
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVSTORE_H
#  include "AliMUONVStore.h"
#endif

class AliMUONLocalTrigger;
class AliMUONGlobalTrigger;
class AliMUONRegionalTrigger;

class AliMUONVTriggerStore : public AliMUONVStore
{
public:
  AliMUONVTriggerStore();
  virtual ~AliMUONVTriggerStore();

  /// Add an object, if of the right type
  virtual Bool_t Add(TObject* object);
  
  /// Add local trigger
  virtual void Add(const AliMUONLocalTrigger& localTrigger) = 0;
  /// Set global trigger
  virtual void SetGlobal(const AliMUONGlobalTrigger& globalTrigger) = 0;
  /// Add regional trigger
  virtual void Add(const AliMUONRegionalTrigger& regionalTrigger) = 0;
  
  using AliMUONVStore::Create;
  
  /// Create a store from the tree (if possible).
  static AliMUONVTriggerStore* Create(TTree& tree);
  
  /// Create iterator (on local card)
  virtual TIterator* CreateIterator() const;
  
  /// Create iterator on local trigger
  virtual TIterator* CreateLocalIterator() const = 0;
  /// Create iterator on regional trigger
  virtual TIterator* CreateRegionalIterator() const = 0;
  
  /// Return global trigger
  virtual AliMUONGlobalTrigger* Global() const = 0;
  
  /// Find a local trigger by the board number (not an index, it is a number really)
  virtual AliMUONLocalTrigger* FindLocal(Int_t boardNumber) const = 0;
  
  /// Find a regional trigger by the board number (not an index, it is a number really)
  virtual AliMUONRegionalTrigger* FindRegional(Int_t boardNumber) const = 0;
  
  using AliMUONVStore::Print;
  
  virtual void Print(Option_t* wildcard="") const { return Print(wildcard,""); }

  ClassDef(AliMUONVTriggerStore,1) // Base class of a trigger store
};

#endif
