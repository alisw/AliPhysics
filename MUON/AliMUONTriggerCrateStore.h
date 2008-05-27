#ifndef ALIMUONTRIGGERCRATESTORE_H
#define ALIMUONTRIGGERCRATESTORE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup trigger
/// \class AliMUONTriggerCrateStore
/// \brief A container for AliMUONTriggerCrate objects.
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#include "TString.h"

class AliMUONLocalTriggerBoard;
class AliMUONTriggerCrate;
class AliMpExMap;
class TIterator;
class AliMUONCalibrationData;

class AliMUONTriggerCrateStore : public TObject
{
public:
  AliMUONTriggerCrateStore();
  virtual ~AliMUONTriggerCrateStore();
  
  Int_t NumberOfCrates() const;

  AliMUONTriggerCrate* Crate(const char* crateName) const;
  AliMUONTriggerCrate* Crate(Int_t ddl, Int_t reg) const;

  Int_t NumberOfLocalBoards() const;

  AliMUONLocalTriggerBoard* LocalBoard(Int_t boardNumber) const;
  
  void ReadFromFile(AliMUONCalibrationData* calibData);
  TIterator* CreateCrateIterator() const;
  
  TIterator* CreateLocalBoardIterator() const;

protected:
  /// Not implemented
  AliMUONTriggerCrateStore(const AliMUONTriggerCrateStore& rhs);
  /// Not implemented
  AliMUONTriggerCrateStore& operator = (const AliMUONTriggerCrateStore& rhs);

private:
  void AddCrate(const char* crateName); 
  
private:
  AliMpExMap* fCrates; ///< list of crates
  AliMpExMap* fLocalBoards; ///< local boards (indexed by their number)

  TString GetCrateName(Int_t ddl, Int_t reg) const;

  ClassDef(AliMUONTriggerCrateStore,2) // Reader for CRATE.TXT file
};

#endif
