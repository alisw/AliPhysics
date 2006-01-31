/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUON2DStore
/// \brief Generic container indexed by a pair of integers
/// 
/// \author Laurent Aphecetche

#ifndef AliMUONV2DSTORE_H
#define AliMUONV2DSTORE_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONV2DStore : public TObject
{
public:
  virtual ~AliMUONV2DStore();
  
  virtual TObject* Get(Int_t i, Int_t j) const = 0;
  virtual Bool_t Set(Int_t i, Int_t j, TObject*, Bool_t replace) = 0;
  virtual Bool_t IsOwner() const = 0;
  
private:  
  ClassDef(AliMUONV2DStore,0) // 
};

#endif
