/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUON1DStore
/// \brief Generic container indexed by a detection element ID
/// 
/// \author Laurent Aphecetche

#ifndef AliMUONV1DSTORE_H
#define AliMUONV1DSTORE_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONV1DStore : public TObject
{
public:
  
  virtual ~AliMUONV1DStore();
  
  virtual TObject* Get(Int_t detElemId) const = 0;

  virtual Bool_t Set(Int_t detElemId, TObject*, Bool_t replace) = 0;
  
  virtual Bool_t IsOwner() const = 0;
  
private:
    
  ClassDef(AliMUONV1DStore,0) // 
};

#endif
