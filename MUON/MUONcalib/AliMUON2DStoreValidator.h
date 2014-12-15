#ifndef ALIMUON2DSTOREVALIDATOR_H
#define ALIMUON2DSTOREVALIDATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUON2DStoreValidator
/// \brief Determine which channels, manus, DEs, stations are missing
/// from a 2DStore.
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVStore;
class TList;
class TObjArray;
class AliMUONCheckItem;
class AliMUONVCalibParam;

class AliMUON2DStoreValidator : public TObject
{
public:
  AliMUON2DStoreValidator();
  virtual ~AliMUON2DStoreValidator();
  
  TObjArray* Validate(const AliMUONVStore& store, Float_t invalidFloatValue, AliMUONVStore* config=0x0);

  TObjArray* Validate(const AliMUONVStore& store, AliMUONVStore* config=0x0);
  
  TObjArray* Validate(const AliMUONVStore& store, 
                      Bool_t (*check)(const AliMUONVCalibParam&,Int_t),
                      AliMUONVStore* config=0x0);

  /// Return statuses
  AliMUONVStore* GetStatus() const { return fStatus; }
  
  /// Reports what is missing, trying to be as concise as possible.
  void Report(TList& lines) const;

  static void Report(TList& lines, const TObjArray& chambers);

private:
  /// Not implemented  
  AliMUON2DStoreValidator(const AliMUON2DStoreValidator&);
  /// Not implemented  
  AliMUON2DStoreValidator& operator=(const AliMUON2DStoreValidator&);

  void AddMissingChannel(Int_t detElemId, Int_t manuId, Int_t manuChannel);

  void AddMissingManu(Int_t detElemId, Int_t manuId);

  AliMUONCheckItem* GetChamber(Int_t chamberID);
  AliMUONCheckItem* GetDE(Int_t detElemId);
  AliMUONCheckItem* GetManu(Int_t detElemId, Int_t manuId);
  
  static void ReportChamber(TList& list, const AliMUONCheckItem& chamber);
  static void ReportDE(TList& list, const AliMUONCheckItem& de);
  static void ReportManu(TList& list, const AliMUONCheckItem& manu);
  
private:
  TObjArray* fChambers; //!< Array of AliMUONCheckItem.
  AliMUONVStore* fStatus; //!< Statuses
  
  ClassDef(AliMUON2DStoreValidator,3) // Validator of 2DStore
};

#endif
