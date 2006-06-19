#ifndef ALIMUONTRIGGERCRATESTORE_H
#define ALIMUONTRIGGERCRATESTORE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONTriggerCrateStore
/// \brief A container for AliMUONTriggerCrate objects.
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONLocalTriggerBoard;
class AliMUONTriggerCrate;
class AliMpExMap;
class TExMapIter;

class AliMUONTriggerCrateStore : public TObject
{
public:
  AliMUONTriggerCrateStore();
  virtual ~AliMUONTriggerCrateStore();
  
  Int_t NumberOfCrates() const;
  void FirstCrate();
  AliMUONTriggerCrate* NextCrate();
  AliMUONTriggerCrate* Crate(const char* crateName) const;
  
  Int_t NumberOfLocalBoards() const;
  void FirstLocalBoard();
  AliMUONLocalTriggerBoard* NextLocalBoard();
  AliMUONLocalTriggerBoard* LocalBoard(Int_t boardNumber) const;
  
  void ReadFromFile(const char* crateFile =
                    "$ALICE_ROOT/MUON/mapping/data/stationTrigger/crate.dat");

protected:
  AliMUONTriggerCrateStore(const AliMUONTriggerCrateStore& rhs);
  AliMUONTriggerCrateStore& operator = (const AliMUONTriggerCrateStore& rhs);

private:
  void AddCrate(const char* crateName); 
  
private:
  AliMpExMap* fCrates; ///< list of crates
  AliMpExMap* fLocalBoards; ///< local boards (indexed by their number)
  TExMapIter* fCrateIterator; //!< iterator for the crate map above
  TExMapIter* fLBIterator; //!< iterator for boards (through crates)
  AliMUONTriggerCrate* fCurrentCrate; //!< used for iterating on local board
  Int_t fCurrentLocalBoard; //!< used for iterating on local board

  ClassDef(AliMUONTriggerCrateStore,1) // Reader for CRATE.TXT file
};

#endif
