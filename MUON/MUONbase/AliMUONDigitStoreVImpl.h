#ifndef ALIMUONDIGITSTOREVIMPL_H
#define ALIMUONDIGITSTOREVIMPL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONDigitStoreVImpl
/// \brief Base implementation of VDigitStore
/// 
// author Laurent Aphecetche

#ifndef ALIMUONVDIGITSTORE_H
#  include "AliMUONVDigitStore.h"
#endif

class TClonesArray;
class AliMUON2DMap;

class AliMUONDigitStoreVImpl : public AliMUONVDigitStore
{
  friend class AliMUONDigitStoreVImplIterator;
  
public:
  AliMUONDigitStoreVImpl(const char* concreteClassName);
  AliMUONDigitStoreVImpl(const AliMUONDigitStoreVImpl& rhs);
  AliMUONDigitStoreVImpl& operator=(const AliMUONDigitStoreVImpl& rhs);
  virtual ~AliMUONDigitStoreVImpl();
  
  /// Whether we can be connected to a TTree
  virtual Bool_t CanConnect() const { return kTRUE; }
  
  virtual Bool_t Connect(TTree& tree, Bool_t alone = kTRUE) const;
  
  virtual void Clear(Option_t* opt="");
  
  virtual AliMUONVDigit* CreateDigit(Int_t detElemId, Int_t manuId,
                                     Int_t manuChannel, Int_t cathode) const = 0;
  
  using AliMUONVDigitStore::Add;
  
  virtual AliMUONVDigit* Add(const AliMUONVDigit& digit, EReplacePolicy replace);
  
  virtual TIterator* CreateIterator() const;
  
  virtual TIterator* CreateIterator(Int_t firstDetElemId, 
                                    Int_t lastDetElemId,
                                    Int_t cathode=2) const;
  
  virtual TIterator* CreateTrackerIterator() const;
  
  virtual TIterator* CreateTriggerIterator() const;
  
  using AliMUONVStore::FindObject;

  virtual AliMUONVDigit* FindObject(UInt_t uniqueID) const;

  virtual AliMUONVDigit* FindObject(Int_t detElemId, Int_t manuId, 
                                    Int_t manuChannel, Int_t cathode) const;
  
  using AliMUONVDigitStore::GetSize;
  
  virtual Int_t GetSize() const;
  
  virtual AliMUONVDigit* Remove(AliMUONVDigit& digit);

protected:
  /// Add concrete digit
  virtual AliMUONVDigit* AddConcreteDigit(TClonesArray& a, 
                                          const AliMUONVDigit& digit,
                                          Int_t index) = 0;
  
private:
  
  AliMUONVDigit* Find(const AliMUONVDigit& digit) const;
  
  void UpdateIndex(const AliMUONVDigit& digit, Int_t index);
  
  Int_t FindIndex(const AliMUONVDigit& digit) const;
  Int_t FindIndex(Int_t detElemId, Int_t internalManuId, Int_t manuChannel) const;
    
  void ReIndex();
  void ClearIndex();

private:
  TClonesArray* fDigits; ///< collection of digits
  AliMUON2DMap* fMap; //!<! index map for fast digit retrieval
  Bool_t fIndexed; //!<! whether our internal indices fDEs and fManus are uptodate
  
  ClassDef(AliMUONDigitStoreVImpl,1) // Implementation of AliMUONVDigitStore
};

#endif
