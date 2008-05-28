#ifndef ALIMUONDIGITSTOREV1ITERATOR_H
#define ALIMUONDIGITSTOREV1ITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONDigitStoreV1Iterator
///
/// \brief Implementation of TIterator for AliMUONDigitStoreV1
/// 
// Author Laurent Aphecetche

#ifndef ALIMUONTOTCASTOREITERATOR_H
#  include "AliMUONTOTCAStoreIterator.h"
#endif

class AliMUONDigitStoreV1Iterator : public AliMUONTOTCAStoreIterator
{
public:
  AliMUONDigitStoreV1Iterator(const AliMUONDigitStoreV1Iterator& rhs);
  AliMUONDigitStoreV1Iterator& operator=(const TIterator& rhs);
  AliMUONDigitStoreV1Iterator& operator=(const AliMUONDigitStoreV1Iterator& rhs);
  AliMUONDigitStoreV1Iterator(TObjArray* a,
                              Int_t firstDetElemId,
                              Int_t lastDetElemId,
                              Int_t cathode=2);
  
  virtual ~AliMUONDigitStoreV1Iterator();
  
  virtual TObject* Next();

  virtual const TCollection* GetCollection() const;
  
private:
  TObjArray* fArray; ///< array we iterate upon
  Int_t fFirstDetElemId; ///< first detection element to iterate upon
  Int_t fLastDetElemId; ///< last detection element to iterate upon
  Int_t fCathode; ///< cathode to iterate upon
  
  ClassDef(AliMUONDigitStoreV1Iterator,1) // Implementation of TIterator
};

#endif
