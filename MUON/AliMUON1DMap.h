/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup
/// \class AliMUONDeMap
/// \brief Basic implementation of AliMUONV1DStore container using
/// AliMpExMap internally
///
/// Should be revised for better preformance ...
/// 
/// \author Laurent Aphecetche

#ifndef AliMUON1DMAP_H
#define AliMUON1DMAP_H

class AliMpExMap;

#include "AliMUONV1DStore.h"

class AliMUON1DMap : public AliMUONV1DStore
{
public:
  AliMUON1DMap();
  virtual ~AliMUON1DMap();

  virtual TObject* Get(Int_t detElemId) const;
  virtual Bool_t Set(Int_t detElemId, TObject* object, Bool_t replace=kTRUE);
  virtual Bool_t IsOwner() const;
  
private:
  AliMpExMap* fMap;
  
  ClassDef(AliMUON1DMap,1) // 
};

#endif
