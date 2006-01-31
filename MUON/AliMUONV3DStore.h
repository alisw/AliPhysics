/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONV3DStore
/// \brief Generic containers indexed by a triplet  
/// (typically detElemId,manuId,manuChannel)
/// 
/// \author Laurent Aphecetche

#ifndef ALIMUONV3DSTORE_H
#define ALIMUONV3DSTORE_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONV3DStore : public TObject
{
public:
  virtual ~AliMUONV3DStore();

  virtual TObject* Get(Int_t i, Int_t j, Int_t k) const = 0;
  virtual Bool_t Set(Int_t i, Int_t j, Int_t k, TObject*, Bool_t replace) = 0;
  virtual Bool_t IsOwner() const = 0;
  
private:
  ClassDef(AliMUONV3DStore,0) // 
};

#endif
