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
// $MpId: AliMpDEIterator.cxx,v 1.6 2006/05/24 13:58:34 ivana Exp $
// Category: management

// ------------------------
// Class AliMpDEIterator
// ------------------------
// The iterator over valid detection elements
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMpDEIterator.h"
#include "AliMpDEStore.h"
#include "AliMpDetElement.h"
#include "AliMpDEManager.h"
#include "AliMpFiles.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TSystem.h>

/// \cond CLASSIMP
ClassImp(AliMpDEIterator)
/// \endcond

//______________________________________________________________________________
AliMpDEIterator::AliMpDEIterator()
    : TObject(),
      fDEStore(AliMpDEStore::Instance()),
      fIndex(-1),
      fChamberId(-1)
{  
/// Standard and default constructor
}

//______________________________________________________________________________
AliMpDEIterator::AliMpDEIterator(const AliMpDEIterator& rhs)
 : TObject(rhs),
   fDEStore(rhs.fDEStore),
   fIndex(rhs.fIndex),
   fChamberId(rhs.fChamberId)
{
/// Copy constructor
}

//______________________________________________________________________________

AliMpDEIterator::~AliMpDEIterator()
{
/// Destructor
}

//______________________________________________________________________________
AliMpDEIterator&  AliMpDEIterator::operator=(const AliMpDEIterator& rhs)
{
/// Assignement operator

  // check assignment to self
  if (this == &rhs) return *this;

  // base class assignment
  TObject::operator=(rhs);

  fDEStore = rhs.fDEStore;
  fIndex = rhs.fIndex;
  fChamberId = rhs.fChamberId;

  return *this;
} 

//
// private methods
//

//______________________________________________________________________________
AliMpDetElement*  AliMpDEIterator::GetDetElement(Int_t index) const
{
/// Return the detection element from the map via index

  return static_cast<AliMpDetElement*>(fDEStore->fDetElements.GetObject(index));
}

//
// public methods
//

//______________________________________________________________________________
void AliMpDEIterator::First()
{
/// Set iterator to the first DE Id defined 

  fIndex = 0;
  fChamberId = -1;
}  

//______________________________________________________________________________
void AliMpDEIterator::First(Int_t chamberId)
{
/// Reset the iterator, so that it points to the first DE
 
  fChamberId = -1;
  fIndex = -1;  
  if ( ! AliMpDEManager::IsValidChamberId(chamberId) ) {
    AliErrorStream() << "Invalid chamber Id " << chamberId << endl;
    return;
  }    

  Int_t i=0;
  while ( i < fDEStore->fDetElements.GetSize() && fChamberId < 0 ) {
    Int_t detElemId = GetDetElement(i)->GetId();
    if ( AliMpDEManager::GetChamberId(detElemId) == chamberId ) {
      fChamberId = chamberId;
      fIndex = i;
    } 
    i++; 
  }

  if ( fChamberId < 0 ) {
    AliErrorStream() 
      << "No DEs of Chamber Id " << chamberId << " found" << endl;
    return;
  }    

}

//______________________________________________________________________________
void AliMpDEIterator::Next()
{
/// Increment iterator to next DE

  fIndex++;

  // Invalidate if at the end
  if ( ( fIndex == fDEStore->fDetElements.GetSize() ) ||
       ( fChamberId >= 0 &&    
         AliMpDEManager::GetChamberId(CurrentDEId()) != fChamberId ) ) {
    fIndex = -1;
  }   
}

//______________________________________________________________________________
Bool_t AliMpDEIterator::IsDone() const
{
/// Is the iterator in the end?

  return ( fIndex < 0 );
}   

//______________________________________________________________________________
AliMpDetElement* AliMpDEIterator::CurrentDE() const
{
/// Current DE Id

  if ( ! IsDone() )
    return GetDetElement(fIndex);
  else {   
    AliErrorStream()
      << "Not in valid position - returning invalid DE." << endl;
    return 0;
  }  
}
    
//______________________________________________________________________________
Int_t AliMpDEIterator::CurrentDEId() const
{
/// Current DE Id

  if ( ! IsDone() )
    return GetDetElement(fIndex)->GetId();
  else {   
    AliErrorStream()
      << "Not in valid position - returning invalid DE." << endl;
    return 0;
  }  
}
    
