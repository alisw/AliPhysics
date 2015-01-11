/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMUON1DArray
/// \brief Implementation of AliMUONVStore
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUON1DARRAY_H
#define ALIMUON1DARRAY_H

#ifndef ALIMUONVSTORE_H
#  include "AliMUONVStore.h"
#endif

class TObjArray;

class AliMUON1DArray : public AliMUONVStore
{
public:
  AliMUON1DArray(Int_t theSize=0);
  AliMUON1DArray(const AliMUON1DArray& other);
  AliMUON1DArray& operator=(const AliMUON1DArray& other);
  
  virtual ~AliMUON1DArray();
  
  virtual AliMUON1DArray* Create() const;
  
  /// Add an object. Object must have a valid UniqueID, which is
  /// used as the index of the array.
  virtual Bool_t Add(TObject* object);

  virtual Bool_t CanConnect() const { return kFALSE; }
  
  virtual void Clear(Option_t* opt="");

  virtual TIterator* CreateIterator() const;
  
  using AliMUONVStore::FindObject;
  
  /// Return the object stored with id.
  virtual TObject* FindObject(UInt_t identifier) const;
    
  using AliMUONVStore::GetSize;
  
  virtual Int_t GetSize() const;
  
private:
   void CopyTo(AliMUON1DArray& to) const;
  Bool_t Set(Int_t i, TObject* object, Bool_t replace);
  
private:  
    
    TObjArray* fArray; ///< Internal array
  
    ClassDef(AliMUON1DArray,1) // Implementation of AliMUONVStore
};

#endif
