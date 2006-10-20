/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONV2DStore
/// \brief Generic container indexed by a pair of integers.
/// 
//  Author Laurent Aphecetche

#ifndef AliMUONV2DSTORE_H
#define AliMUONV2DSTORE_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVDataIterator;

class AliMUONV2DStore : public TObject
{
public:
  AliMUONV2DStore();
  virtual ~AliMUONV2DStore();
  
  /// Return an empty copy of self.
  virtual AliMUONV2DStore* CloneEmpty() const { return 0x0; }
  
  virtual AliMUONVDataIterator* Iterator() const { return 0x0; }
  
  /// Return the object stored at (i,j).
  virtual TObject* Get(Int_t i, Int_t j) const = 0;
  
  /** Set the object stored at (i,j).
    if replace=false and there's already an object there, returns kFALSE
    */
  virtual Bool_t Set(Int_t i, Int_t j, TObject*, Bool_t replace) = 0;
  
  /// Whether or not this container is the owner of its contents.
  virtual Bool_t IsOwner() const = 0;
  
private:  
  ClassDef(AliMUONV2DStore,0) // Generic container indexed by a pair of integers
};

#endif
