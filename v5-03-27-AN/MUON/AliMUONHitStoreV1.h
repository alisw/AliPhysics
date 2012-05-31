#ifndef ALIMUONHITSTOREV1_H
#define ALIMUONHITSTOREV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup sim
/// \class AliMUONHitStoreV1
/// \brief Implementation of AliMUONVHitStore
/// 
// Author Laurent Aphecetche

#ifndef ALIMUONVHITSTORE_H
#  include "AliMUONVHitStore.h"
#endif

class TClonesArray;

class AliMUONHitStoreV1 : public AliMUONVHitStore
{
public:
  AliMUONHitStoreV1();
  AliMUONHitStoreV1(TRootIOCtor* /*dummy*/);
  virtual ~AliMUONHitStoreV1();
  
  using AliMUONVHitStore::Add;

  virtual void Add(const AliMUONHit& hit);

  /// Whether the Connect(TTree&) method is implemented
  virtual Bool_t CanConnect() const { return kTRUE; }

  virtual void Clear(Option_t* opt="");
  
  virtual Bool_t Connect(TTree& tree, Bool_t alone=kTRUE) const;
  
  virtual AliMUONHitStoreV1* Create() const { return new AliMUONHitStoreV1; }
  
  virtual TIterator* CreateIterator() const;
  
  virtual TCollection* Collection();
  
  using AliMUONVHitStore::GetSize;
  
  virtual Int_t GetSize() const;
  
private:
  /// Not implemented
  AliMUONHitStoreV1(const AliMUONHitStoreV1&);
  /// Not implemented
  AliMUONHitStoreV1& operator=(const AliMUONHitStoreV1&);
  /// Return the address of array of hits
  TClonesArray** HitsPtr() const { return const_cast<TClonesArray**>(&fHits); }
  /// Return the array of hits
  TClonesArray* Hits() const { return fHits; }

private:
    TClonesArray* fHits; ///< array of hits
  
  ClassDef(AliMUONHitStoreV1,1) // Implementation of AliMUONVHitStore
};

#endif
