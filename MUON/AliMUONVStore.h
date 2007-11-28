#ifndef ALIMUONVSTORE_H
#define ALIMUONVSTORE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMUONVStore
/// \brief Base class for MUON data stores.
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TIterator;
class TTree;

class AliMUONVStore : public TObject
{
public:
  AliMUONVStore();
  virtual ~AliMUONVStore();
  
  /// Add an object to the store
  virtual Bool_t Add(TObject* object) = 0;
  
  /// Clear ourselves (i.e. Reset)
  virtual void Clear(Option_t* opt="") = 0;
  
  /// Create an empty copy of this
  virtual AliMUONVStore* Create() const = 0;

  /// Create a store from a TTree
  static AliMUONVStore* Create(TTree& tree, const char* what);

  /// Return an iterator to loop over the whole store
  virtual TIterator* CreateIterator() const = 0;
  
  /// Whether the Connect(TTree&) method is implemented
  virtual Bool_t CanConnect() const = 0;
  
  /// Connect us to a TTree (only valid if CanConnect()==kTRUE)
  virtual Bool_t Connect(TTree& tree, Bool_t alone=kTRUE) const;

  /// Find an object by name
  virtual TObject* FindObject(const char* name) const;
  
  /// Find an object
  virtual TObject* FindObject(const TObject* object) const;

  /// Find an object using a single id
	virtual TObject* FindObject(UInt_t uniqueID) const;
  
  /// Find an object using 2 ids
  virtual TObject* FindObject(Int_t i, Int_t j) const;
  
  /// The number of objects stored
  virtual Int_t GetSize() const = 0;

  /// The number of objects stored for firstid=i. Not implemented by default.
  virtual Int_t GetSize(Int_t i) const;

  /// Whether we are empty or not
  virtual Bool_t IsEmpty() const { return GetSize() == 0; }
  
  /// Print all objects whose name matches wildcard
  virtual void Print(Option_t* wildcard="") const;

  /// Print, with option, all objects whose name matches wildcard
  virtual void Print(Option_t* wildcard, Option_t* opt) const;
  
  ClassDef(AliMUONVStore,1) // Base class for a MUON data store
};

#endif
