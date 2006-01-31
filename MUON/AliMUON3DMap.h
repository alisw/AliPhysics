/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUON3DMap
/// \brief Basic implementation of AliMUONV3DStore container using
/// AliMpExMap internally
///
/// Should be revised for better preformance ...
/// 
/// \author Laurent Aphecetche

#ifndef ALIMUON3DMAP_H
#define ALIMUON3DMAP_H

#ifndef ALIMUONV3DSTORE_H
#  include "AliMUONV3DStore.h"
#endif

class AliMUONV1DStore;

class AliMUON3DMap : public AliMUONV3DStore
{
public:
  AliMUON3DMap();
  virtual ~AliMUON3DMap();
  
  virtual TObject* Get(Int_t i, Int_t j, Int_t k) const;
  virtual Bool_t Set(Int_t i, Int_t j, Int_t k, TObject*, Bool_t replace);
  virtual Bool_t IsOwner() const;
  virtual void Print(Option_t* opt="") const;
  
private:  
    AliMUONV1DStore* fStore;  
  ClassDef(AliMUON3DMap,1) // 
};

#endif
