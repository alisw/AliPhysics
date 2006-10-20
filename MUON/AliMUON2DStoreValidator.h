#ifndef ALIMUON2DSTOREVALIDATOR_H
#define ALIMUON2DSTOREVALIDATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUON2DStoreValidator
/// \brief Determine which channels, manus, DEs, stations are missing
/// from a 2DStore.
/// 
/// \author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONV2DStore;
class TList;
class TObjArray;
class AliMUONCheckItem;

class AliMUON2DStoreValidator : public TObject
{
public:
  AliMUON2DStoreValidator();
  virtual ~AliMUON2DStoreValidator();
  
  TObjArray* Validate(const AliMUONV2DStore& store, Float_t invalidFloatValue);
  void Report() const { Report(*fChambers); }

  static void Report(const TObjArray& chambers);

private:
    
  AliMUON2DStoreValidator(const AliMUON2DStoreValidator&);
  AliMUON2DStoreValidator& operator=(const AliMUON2DStoreValidator&);

  void AddMissingChannel(Int_t detElemId, Int_t manuId, Int_t manuChannel);

  void AddMissingManu(Int_t detElemId, Int_t manuId);

  AliMUONCheckItem* GetChamber(Int_t chamberID);
  AliMUONCheckItem* GetDE(Int_t detElemId);
  AliMUONCheckItem* GetManu(Int_t detElemId, Int_t manuId);
  
  static void ReportChamber(AliMUONCheckItem& chamber);
  static void ReportDE(AliMUONCheckItem& de);
  static void ReportManu(AliMUONCheckItem& manu);
  
private:
  TList* fManuList; //! List of (DE,manuID) pairs.
  TObjArray* fChambers; //! Array of AliMUONCheckItem.
  
  ClassDef(AliMUON2DStoreValidator,1) // Validator of 2DStore
};

#endif
