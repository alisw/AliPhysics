/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONClusterStoreV2Iterator
///
/// Implementation of TIterator for AliMUONClusterStoreV2
///
/// \author Philippe Pillot, Subatech
///
//-----------------------------------------------------------------------------

#include "AliMUONClusterStoreV2Iterator.h"

#include "AliMUONClusterStoreV2.h"

#include "AliMpExMap.h"

#include <TExMap.h>

/// \cond CLASSIMP
ClassImp(AliMUONClusterStoreV2Iterator)
/// \endcond

//_____________________________________________________________________________
AliMUONClusterStoreV2Iterator::AliMUONClusterStoreV2Iterator(const AliMUONClusterStoreV2* store,
                                                             Int_t firstChamberId, Int_t lastChamberId)
: TIterator(),
  fStore(store),
  fFirstChamberId(firstChamberId),
  fLastChamberId(lastChamberId),
  fCurrentChamberId(-1),
  fChamberIterator(0x0)
{
  /// Constructor for partial iteration
  if (fFirstChamberId > fLastChamberId) {
    fLastChamberId = fFirstChamberId;
    fFirstChamberId = lastChamberId;
  }
  Reset();
}

//_____________________________________________________________________________
AliMUONClusterStoreV2Iterator::AliMUONClusterStoreV2Iterator(const AliMUONClusterStoreV2Iterator& iter)
: TIterator(iter),
  fStore(iter.fStore),
  fFirstChamberId(iter.fFirstChamberId),
  fLastChamberId(iter.fLastChamberId),
  fCurrentChamberId(iter.fCurrentChamberId),
  fChamberIterator(0x0)
{
  /// Copy constructor
  if (iter.fChamberIterator) fChamberIterator = new TExMapIter(*(iter.fChamberIterator));
}

//_____________________________________________________________________________
AliMUONClusterStoreV2Iterator& AliMUONClusterStoreV2Iterator::operator=(const AliMUONClusterStoreV2Iterator& iter)
{
  /// Assignment operator
  if ( this != &iter ) {
    fFirstChamberId = iter.fFirstChamberId;
    fLastChamberId = iter.fLastChamberId;
    fCurrentChamberId = iter.fCurrentChamberId;
    delete fChamberIterator;
    fChamberIterator = 0x0;
    if (iter.fChamberIterator) fChamberIterator = new TExMapIter(*(iter.fChamberIterator));
  }
  return *this;
}

//_____________________________________________________________________________
TIterator& AliMUONClusterStoreV2Iterator::operator=(const TIterator& iter)
{
  /// Overriden operator= (imposed by Root's definition of TIterator::operator= ?)
  if ( this != &iter && iter.IsA() == AliMUONClusterStoreV2Iterator::Class() ) {
    (*this) = static_cast<const AliMUONClusterStoreV2Iterator&>(iter);
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONClusterStoreV2Iterator::~AliMUONClusterStoreV2Iterator()
{
  /// Destructor
  delete fChamberIterator;
}

//_____________________________________________________________________________
TObject* AliMUONClusterStoreV2Iterator::NextInCurrentChamber() const
{
  /// Return the value corresponding to theKey in iterator iter
  Long_t key, value;
  if (fChamberIterator->Next(key,value)) return reinterpret_cast<TObject*>(value);
  else return 0x0;
}

//_____________________________________________________________________________
TObject* AliMUONClusterStoreV2Iterator::Next()
{
  /// Return next cluster in store
  TObject* o = NextInCurrentChamber();
  
  while (!o) {
    // fChamberIterator exhausted, try to get the next ones
    if (fCurrentChamberId == fLastChamberId) return 0x0; // we reached the end
    
    fCurrentChamberId++;
    delete fChamberIterator;
    fChamberIterator = new TExMapIter(static_cast<AliMpExMap*>(fStore->fMap->UncheckedAt(fCurrentChamberId))->GetIterator());
    
    o = NextInCurrentChamber();
  }
  
  return o;
}

//_____________________________________________________________________________
void AliMUONClusterStoreV2Iterator::Reset()
{
  /// Reset the iterator
  fCurrentChamberId = fFirstChamberId;
  delete fChamberIterator;
  fChamberIterator = new TExMapIter(static_cast<AliMpExMap*>(fStore->fMap->UncheckedAt(fCurrentChamberId))->GetIterator());
}
