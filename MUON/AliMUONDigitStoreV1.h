#ifndef ALIMUONDIGITSTOREV1_H
#define ALIMUONDIGITSTOREV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONDigitStoreV1
/// \brief (Legacy) implementation of AliMUONVDigitStore
/// 
// Author Laurent Aphecetche

#ifndef ALIMUONVDIGITSTORE_H
#  include "AliMUONVDigitStore.h"
#endif

class TObjArray;
class TClonesArray;

class AliMUONDigitStoreV1 : public AliMUONVDigitStore
{
public:
  AliMUONDigitStoreV1();
  AliMUONDigitStoreV1(const AliMUONDigitStoreV1& rhs);
  AliMUONDigitStoreV1& operator=(const AliMUONDigitStoreV1& rhs);  
  virtual ~AliMUONDigitStoreV1();
  
  virtual void Clear(Option_t* opt="");

  /// Whether the Connect(TTree&) method is implemented
  virtual Bool_t CanConnect() const { return kTRUE; }
  
  // Below are our specific methods
    
  virtual AliMUONVDigitStore* Create() const { return new AliMUONDigitStoreV1; }
  
  virtual Bool_t Connect(TTree& tree, Bool_t alone=kTRUE) const;

  virtual AliMUONVDigit* CreateDigit(Int_t detElemId, Int_t manuId,
                                     Int_t manuChannel, Int_t cathode) const;
  
  using AliMUONVDigitStore::Add;
  
  virtual AliMUONVDigit* Add(const AliMUONVDigit& digit, EReplacePolicy replace=kDeny);

  virtual TIterator* CreateIterator() const;
  
  virtual TIterator* CreateTrackerIterator() const;
  
  virtual TIterator* CreateTriggerIterator() const;

  virtual TIterator* CreateIterator(Int_t firstDetElemId, 
                                    Int_t lastDetElemId,
                                    Int_t cathode=2) const;
    
  using AliMUONVDigitStore::FindObject;
  
  virtual AliMUONVDigit* FindObject(Int_t detElemId, Int_t manuId, 
                                    Int_t manuChannel, Int_t cathode) const;

  using AliMUONVDigitStore::GetSize;
  
  virtual Int_t GetSize() const;
  
  virtual AliMUONVDigit* Remove(AliMUONVDigit& digit);

  Bool_t HasMCInformation() const;
  
private:

  TObject** ChamberDigitsPtr(Int_t chamberId) const;
  
  TClonesArray* ChamberDigits(Int_t chamberId);
  const TClonesArray* ChamberDigits(Int_t chamberId) const;
  
  AliMUONVDigit* Find(const AliMUONVDigit& digit, Int_t& index) const;
  
  AliMUONVDigit* FindIndex(Int_t detElemId, Int_t manuId, 
                           Int_t manuChannel, Int_t cathode, Int_t& index) const;
  
private:
  TObjArray* fDigits; ///< array of tclonesarray
  TClonesArray* fChamberDigits; ///< array of digits for one chamber
  
  ClassDef(AliMUONDigitStoreV1,1) // (Legacy) Implementation of AliMUONVDigitStore
};

#endif
