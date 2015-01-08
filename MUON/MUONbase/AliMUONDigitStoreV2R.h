#ifndef ALIMUONDIGITSTOREV2R_H
#define ALIMUONDIGITSTOREV2R_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONDigitStoreV2R
/// \brief Concrete implementation of AliMUONVDigitStore for real digits
/// 
// author Laurent Aphecetche

#ifndef ALIMUONDIGITSTOREVIMPL_H
#  include "AliMUONDigitStoreVImpl.h"
#endif

class AliMUONDigitStoreV2R : public AliMUONDigitStoreVImpl
{
public:
  AliMUONDigitStoreV2R();
  virtual ~AliMUONDigitStoreV2R();
  
  virtual AliMUONVDigitStore* Create() const;

  virtual AliMUONVDigit* CreateDigit(Int_t detElemId, Int_t manuId,
                                     Int_t manuChannel, Int_t cathode) const;

  virtual Bool_t HasMCInformation() const { return kFALSE; }
  
protected:
    
    virtual AliMUONVDigit* AddConcreteDigit(TClonesArray& a, 
                                            const AliMUONVDigit& digit,
                                            Int_t index);

  ClassDef(AliMUONDigitStoreV2R,1) // Implementation of AliMUONVDigitStore
};

#endif
