#ifndef ALIMUONDIGITSTOREVIMPLITERATOR_H
#define ALIMUONDIGITSTOREVIMPLITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONDigitStoreVImplIterator
/// \brief Base implementation of TIterator for AliMUONDigitStoreVImpl
/// 
// Author Laurent Aphecetche

#ifndef ROOT_TIterator
#  include "TIterator.h"
#endif

class AliMUONDigitStoreVImpl;
class AliMUONVCalibParam;

class AliMUONDigitStoreVImplIterator : public TIterator
{
public:
  AliMUONDigitStoreVImplIterator(const AliMUONDigitStoreVImpl* store);
  AliMUONDigitStoreVImplIterator(const AliMUONDigitStoreVImpl* store,
                              Int_t firstDetElemId,
                              Int_t lastDetElemId,
                              Int_t cathode=2);
  
  virtual ~AliMUONDigitStoreVImplIterator();
  
  TObject* Next();
  
  void Reset();
  
  /// Return 0 as we're not dealing with TCollection objects really
  virtual const TCollection* GetCollection() const { return 0x0; }
  
private:
  /// Not implemented
  AliMUONDigitStoreVImplIterator(const AliMUONDigitStoreVImplIterator& rhs);
  /// Not implemented
  AliMUONDigitStoreVImplIterator& operator=(const AliMUONDigitStoreVImplIterator& rhs);
  /// Overriden TIterator virtual operator=
  AliMUONDigitStoreVImplIterator& operator=(const TIterator& rhs);

  const AliMUONDigitStoreVImpl* fkStore; //!<! store to iterate upon
  Int_t fFirstDetElemId; //!<! first de
  Int_t fLastDetElemId; //!<! last de
  Int_t fCathode; //!<! cathode (-1 for both)
  TIterator* fStoreIterator; //!<! helper iterator
  AliMUONVCalibParam* fCurrentCalibParam; //!<! current CalibParam
  Int_t fCurrentCalibParamIndex; //!<! current index in fCurrentCalibParam
  
  ClassDef(AliMUONDigitStoreVImplIterator,1) // Implementation of AliMUONVDataIterator
};

#endif
