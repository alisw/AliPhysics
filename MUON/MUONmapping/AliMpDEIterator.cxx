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

//-----------------------------------------------------------------------------
// Class AliMpDEIterator
// ------------------------
// The iterator over valid detection elements
// Author: Ivana Hrivnacova, IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpDEIterator.h"

#include "AliMpExMapIterator.h"
#include "AliMpDEStore.h"
#include "AliMpDetElement.h"
#include "AliMpDEManager.h"
#include "AliMpFiles.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TSystem.h>

using std::endl;
/// \cond CLASSIMP
ClassImp(AliMpDEIterator)
/// \endcond

//______________________________________________________________________________
AliMpDEIterator::AliMpDEIterator()
: TObject(),
  fCurrentDE(0x0),
  fIterator(AliMpDEStore::Instance()->fDetElements.CreateIterator()),
  fChamberId(-1)
{  
/// Standard and default constructor
}

//______________________________________________________________________________

AliMpDEIterator::~AliMpDEIterator()
{
/// Destructor

  delete fIterator;
}

//
// public methods
//

//______________________________________________________________________________
void AliMpDEIterator::First()
{
/// Set iterator to the first DE Id defined 

  fIterator->Reset();
  fCurrentDE = static_cast<AliMpDetElement*>(fIterator->Next());
  fChamberId = -1;
}  

//______________________________________________________________________________
void AliMpDEIterator::First(Int_t chamberId)
{
/// Reset the iterator, so that it points to the first DE
 
  if ( ! AliMpDEManager::IsValidChamberId(chamberId) ) {
    AliErrorStream() << "Invalid chamber Id " << chamberId << endl;
    fIterator->Reset();
    fChamberId = -1;    
    fCurrentDE = 0x0;
    return;
  }    

  fIterator->Reset();
  fChamberId = -1;
  while ( fChamberId != chamberId ) 
  {
    fCurrentDE = static_cast<AliMpDetElement*>(fIterator->Next());
    if (!fCurrentDE) return;
    fChamberId = AliMpDEManager::GetChamberId(CurrentDEId());
  }
}

//______________________________________________________________________________
void AliMpDEIterator::Next()
{
/// Increment iterator to next DE

  if ( fChamberId < 0 ) 
  {
    fCurrentDE = static_cast<AliMpDetElement*>(fIterator->Next());
  }
  else
  {
    fCurrentDE = static_cast<AliMpDetElement*>(fIterator->Next());
    
    while ( fCurrentDE && (AliMpDEManager::GetChamberId(fCurrentDE->GetId()) != fChamberId) )
    {
      fCurrentDE = static_cast<AliMpDetElement*>(fIterator->Next());
      if (!fCurrentDE) return;
    }
  }
}

//______________________________________________________________________________
Bool_t AliMpDEIterator::IsDone() const
{
/// Is the iterator in the end?

  return ( fCurrentDE == 0x0 );
}   

//______________________________________________________________________________
AliMpDetElement* AliMpDEIterator::CurrentDE() const
{
/// Current DE Id

  return fCurrentDE;
}
    
//______________________________________________________________________________
Int_t 
AliMpDEIterator::CurrentDEId() const
{
/// Current DE Id

  if ( fCurrentDE )
  {
    return fCurrentDE->GetId();
  }
  else {   
    AliErrorStream()
      << "Not in valid position - returning invalid DE." << endl;
    return 0;
  }  
}
    
