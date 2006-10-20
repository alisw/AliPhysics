/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUON2DMap
/// \brief Basic implementation of AliMUONV2DStore container using
/// AliMpExMap internally.
///
//  Author Laurent Aphecetche

#ifndef AliMUON2DMAP_H
#define AliMUON2DMAP_H

#include "AliMUONV2DStore.h"

class AliMpExMap;

class AliMUON2DMap : public AliMUONV2DStore
{
public:
  AliMUON2DMap(Bool_t optimizeForDEManu=kFALSE);  
  virtual ~AliMUON2DMap();

  AliMUONV2DStore* CloneEmpty() const;
  
  /// The returned iterator is owned by the client.
  AliMUONVDataIterator* Iterator() const;
  
  virtual TObject* Get(Int_t i, Int_t j) const;
  virtual Bool_t Set(Int_t i, Int_t j, TObject* object, Bool_t replace);
  /// Whether or not this container is the owner of its contents.
  virtual Bool_t IsOwner() const { return kTRUE; } 

  virtual void Print(Option_t* opt="") const;
  
  AliMUON2DMap(const AliMUON2DMap& other);
  AliMUON2DMap&  operator = (const AliMUON2DMap& other);

private:
  void CopyTo(AliMUON2DMap& destination) const;

private:
  AliMpExMap* fMap; ///< Our internal map (an AliMpExMap of AliMpExMaps)
  Bool_t fOptimizeForDEManu; ///< whether (i,j) pair is supposed to be (DetElemId,ManuId) (allow us to allocate right amount of memory, that's all it does.
  
  ClassDef(AliMUON2DMap,2) // A 2D container
};

#endif
