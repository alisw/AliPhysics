#ifndef ALIMUONTOTCASTOREITERATOR_H
#define ALIMUONTOTCASTOREITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONTOTCAStoreIterator
/// \brief Iterator on a store composed of a TObjArray of TClonesArrays
/// 
//  Author: Laurent Aphecetche

#include "TIterator.h"

class TClonesArray;
class TObjArray;

class AliMUONTOTCAStoreIterator : public TIterator
{
public:
  AliMUONTOTCAStoreIterator(const TObjArray* a, Int_t firstChamberId, Int_t lastChamberId);
  AliMUONTOTCAStoreIterator(const AliMUONTOTCAStoreIterator& rhs);
  AliMUONTOTCAStoreIterator& operator=(const TIterator& rhs);
  AliMUONTOTCAStoreIterator& operator=(const AliMUONTOTCAStoreIterator& rhs);
  virtual ~AliMUONTOTCAStoreIterator();
    
  virtual const TCollection* GetCollection() const;
  
  virtual TObject* Next();
  
  virtual void Reset(); 
  
private:
    void CopyTo(AliMUONTOTCAStoreIterator& destination) const;
  
private:
  const TObjArray* fkData; //!<! Pointer to data accessor
  Int_t fFirstChamberId;      //!<! First chamber to iterate on
  Int_t fLastChamberId;       //!<! Last chamber to iterate on      
  TClonesArray* fCurrentTCA;    //!<! TClonesArray of the current chamber
  Int_t fCurrentTCAIndex;      //!<! Current position within fCurrentTCA array
  Int_t fCurrentChamberId; //!<! current chamber id
  
  ClassDef(AliMUONTOTCAStoreIterator,0) // Iterator on digits
};      

#endif
