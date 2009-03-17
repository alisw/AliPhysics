#ifndef ALIMUONDIGITSTOREV2S_H
#define ALIMUONDIGITSTOREV2S_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONDigitStoreV2S
/// \brief Implementation of AliMUONVDigitStore
/// 
// author Laurent Aphecetche

#ifndef ALIMUONDIGITSTOREVIMPL_H
#  include "AliMUONDigitStoreVImpl.h"
#endif

class AliMUONDigitStoreV2S : public AliMUONDigitStoreVImpl
{
public:
  AliMUONDigitStoreV2S();
  virtual ~AliMUONDigitStoreV2S();
  
  virtual AliMUONVDigitStore* Create() const;

  virtual AliMUONVDigit* CreateDigit(Int_t detElemId, Int_t manuId,
                                     Int_t manuChannel, Int_t cathode) const;

  virtual Bool_t HasMCInformation() const { return kTRUE; }

protected:
    
    virtual AliMUONVDigit* AddConcreteDigit(TClonesArray& a, 
                                            const AliMUONVDigit& digit,
                                            Int_t index);

  ClassDef(AliMUONDigitStoreV2S,1) // Implementation of AliMUONVDigitStore
};

#endif
