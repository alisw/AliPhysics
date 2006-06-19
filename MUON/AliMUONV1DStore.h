/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONV1DStore
/// \brief Generic container indexed by a single integer.
/// 
//  Author Laurent Aphecetche

#ifndef AliMUONV1DSTORE_H
#define AliMUONV1DSTORE_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONV1DStore : public TObject
{
public:
  AliMUONV1DStore();
  virtual ~AliMUONV1DStore();
  
  /// Return the object stored at i.
  virtual TObject* Get(Int_t i) const = 0;
  
  /** Set the object stored at i.
    if replace=false and there's already an object there, returns kFALSE
    */
  virtual Bool_t Set(Int_t i, TObject*, Bool_t replace) = 0;
  
  /// Whether or not this container is the owner of its contents.
  virtual Bool_t IsOwner() const = 0;
  
private:  
  ClassDef(AliMUONV1DStore,0) // Generic container indexed by a single integer
};

#endif
