#ifndef ALIMUONTRIGGERSTOREV1_H
#define ALIMUONTRIGGERSTOREV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONTriggerStoreV1
/// \brief Implementation of AliMUONVTriggerStore
/// 
// Author Laurent Aphecetche

#ifndef ALIMUONVTRIGGERSTORE_H
#  include "AliMUONVTriggerStore.h"
#endif

class TClonesArray;

class AliMUONTriggerStoreV1 : public AliMUONVTriggerStore
{
public:
  AliMUONTriggerStoreV1();
  virtual ~AliMUONTriggerStoreV1();
  
  /// Whether the Connect(TTree&) method is implemented
  virtual AliMUONTriggerStoreV1* Create() const { return new AliMUONTriggerStoreV1; }
  
  virtual Bool_t CanConnect() const { return kTRUE; }
  
  using AliMUONVTriggerStore::Add;
  virtual void Add(const AliMUONLocalTrigger& localTrigger);
  virtual void Add(const AliMUONRegionalTrigger& regionalTrigger);
  
  virtual void SetGlobal(const AliMUONGlobalTrigger& globalTrigger);
  
  virtual TIterator* CreateLocalIterator() const;
  virtual TIterator* CreateRegionalIterator() const;
  
  virtual AliMUONGlobalTrigger* Global() const;

  virtual AliMUONLocalTrigger* FindLocal(Int_t boardNumber) const;
  virtual AliMUONRegionalTrigger* FindRegional(Int_t boardNumber) const;

  using AliMUONVTriggerStore::Print;
  
  virtual void Print(Option_t* wildcard, Option_t* opt) const;

  virtual Bool_t Connect(TTree& tree, Bool_t alone=kTRUE) const;
  
  virtual void Clear(Option_t* opt="");

  using AliMUONVStore::GetSize;
  
  virtual Int_t GetSize() const;
  
private:
  /// Return the address of the array of global trigger information
  TClonesArray** GlobalPtr() const { return const_cast<TClonesArray**>(&fGlobal); }
  /// Return the address of the array of local trigger information
  TClonesArray** LocalPtr() const { return const_cast<TClonesArray**>(&fLocal); }
  /// Return the address of the array of regional trigger information
  TClonesArray** RegionalPtr() const { return const_cast<TClonesArray**>(&fRegional); }
  
  /// Not implemented
  AliMUONTriggerStoreV1(const AliMUONTriggerStoreV1&);
  /// Not implemented
  AliMUONTriggerStoreV1& operator=(const AliMUONTriggerStoreV1&);
  
private:
  TClonesArray* fLocal; ///< internal array of local trigger information
  TClonesArray* fRegional; ///< internal array of regional trigger information
  TClonesArray* fGlobal; ///< internal array of global trigger information
  mutable TClonesArray* fEmptyLocal; //!<! internal array of empty local trigger
  
  ClassDef(AliMUONTriggerStoreV1,1) // Implementation of AliMUONVTriggerStore
};

#endif
