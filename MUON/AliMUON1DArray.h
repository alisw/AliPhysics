/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUON1DArray
/// \brief Implementation of AliMUONV1DStore
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUON1DARRAY_H
#define ALIMUON1DARRAY_H

#ifndef ALIMUONV1DSTORE_H
#  include "AliMUONV1DStore.h"
#endif

class TObjArray;

class AliMUON1DArray : public AliMUONV1DStore
{
public:
  AliMUON1DArray(Int_t theSize=0);
  AliMUON1DArray(const AliMUON1DArray& other);
  AliMUON1DArray& operator=(const AliMUON1DArray& other);
  
  virtual ~AliMUON1DArray();
  
  /// Return the object stored at i.
  virtual TObject* Get(Int_t i) const;
  
  /** Set the object stored at i.
    if replace=false and there's already an object there, returns kFALSE
    */
  virtual Bool_t Set(Int_t i, TObject* object, Bool_t replace);
  
  /// Whether or not this container is the owner of its contents.
  virtual Bool_t IsOwner() const { return kTRUE; }
  
private:
   void CopyTo(AliMUON1DArray& to) const;
  
private:  
    
    TObjArray* fArray; ///< Internal array
  
    ClassDef(AliMUON1DArray,1) // Implementation of AliMUONV1DStore
};

#endif
