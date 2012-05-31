/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMUON2DMap
/// \brief Basic implementation of AliMUONVStore container using
/// AliMpExMap internally.
///
//  Author Laurent Aphecetche

#ifndef ALIMUON2DMAP_H
#define ALIMUON2DMAP_H

#include "AliMUONVStore.h"

class AliMpExMap;

class AliMUON2DMap : public AliMUONVStore
{
public:
  AliMUON2DMap(TRootIOCtor* ioCtor);
  AliMUON2DMap(Bool_t optimizeForDEManu);  
  AliMUON2DMap(const AliMUON2DMap& other);
  AliMUON2DMap&  operator = (const AliMUON2DMap& other);
  virtual ~AliMUON2DMap();

  virtual Bool_t Add(TObject* object);
  
  /// Mandatory methods from TCollection
  virtual void Clear(Option_t* opt="");
  
  /// Whether the Connect(TTree&) method is implemented
  virtual Bool_t CanConnect() const { return kFALSE; }
  
  virtual AliMUONVStore* Create() const;
  
  /// The returned iterator is owned by the client.
  virtual TIterator* CreateIterator() const;

  /// Iterate on part of the store (only for (i,j) where firstI<=i<=lastI
  TIterator* CreateIterator(Int_t firstI, Int_t lastI) const;

  using AliMUONVStore::FindObject;
  
  virtual TObject* FindObject(Int_t i, Int_t j) const;

  virtual TObject* FindObject(UInt_t uid) const;

  /// Whether our internal storage is optimize to store (detection element id, manu id)  
  Bool_t IsOptimizedForDEManu() const { return fOptimizeForDEManu; }
  
  virtual Int_t GetSize() const;

  virtual Int_t GetSize(Int_t i) const;

private:
  /// Not implemented
  AliMUON2DMap();

  Bool_t Set(Int_t i, Int_t j, TObject* object, Bool_t replace);

  static const Int_t fgkOptimalSizeForDEManu; ///< Optimal size for DEManu

  AliMpExMap* fMap; ///< Our internal map (an AliMpExMap of AliMpExMaps)
  Bool_t fOptimizeForDEManu; ///< whether (i,j) pair is supposed to be (DetElemId,ManuId) (allow us to allocate right amount of memory, that's all it does.
  
  ClassDef(AliMUON2DMap,2) // A 2D container
};

#endif
