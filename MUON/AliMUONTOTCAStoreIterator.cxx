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
/// \class AliMUONTOTCAStoreIterator
///
/// An iterator to access TObject stored in a TObjArray of TClonesArray
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONTOTCAStoreIterator.h"

#include "AliLog.h"
#include <TClonesArray.h>
#include <TObjArray.h>

/// \cond CLASSIMP
ClassImp(AliMUONTOTCAStoreIterator)
/// \endcond

//_____________________________________________________________________________
AliMUONTOTCAStoreIterator::AliMUONTOTCAStoreIterator(const TObjArray* data,
                                                     Int_t firstChamberId, 
                                                     Int_t lastChamberId)
: 
TIterator(),
fkData(data),
fFirstChamberId(firstChamberId),
fLastChamberId(lastChamberId),
fCurrentTCA(0x0),
fCurrentTCAIndex(-1),
fCurrentChamberId(-1)
{
  /// Standard constructor
  Reset();
}

//_____________________________________________________________________________
AliMUONTOTCAStoreIterator& 
AliMUONTOTCAStoreIterator::operator=(const TIterator& rhs)
{
  /// Overriden operator= (imposed by Root's declaration of TIterator ?)
  if ( this != &rhs )
  {
    if ( rhs.IsA() != AliMUONTOTCAStoreIterator::Class() )
    {
      AliErrorGeneral("AliMUONTOTCAStoreIterator::operator=","Wrong type");
    }
    else
    {
      const AliMUONTOTCAStoreIterator& rhs1 = 
      static_cast<const AliMUONTOTCAStoreIterator&>(rhs);
      rhs1.CopyTo(*this);
    }
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONTOTCAStoreIterator::AliMUONTOTCAStoreIterator(const AliMUONTOTCAStoreIterator& rhs)
: 
TIterator(rhs),
fkData(0x0),
fFirstChamberId(-1),
fLastChamberId(-1),
fCurrentTCA(0x0),
fCurrentTCAIndex(-1),
fCurrentChamberId(-1)
{
  /// Copy constructor

  rhs.CopyTo(*this);
}

//_____________________________________________________________________________
AliMUONTOTCAStoreIterator::~AliMUONTOTCAStoreIterator()
{
  /// Destructor
}

//_____________________________________________________________________________
AliMUONTOTCAStoreIterator&
AliMUONTOTCAStoreIterator::operator=(const AliMUONTOTCAStoreIterator& rhs)
{
  /// Assignment operator

  rhs.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
void
AliMUONTOTCAStoreIterator::CopyTo(AliMUONTOTCAStoreIterator& destination) const
{
  /// Copy *this to destination
  destination.fkData=fkData;
  destination.fFirstChamberId=fFirstChamberId;
  destination.fLastChamberId=fLastChamberId;
  destination.fCurrentTCAIndex=fCurrentTCAIndex;
  destination.fCurrentChamberId=fCurrentChamberId;
  destination.fCurrentTCA=fCurrentTCA;
}

//_____________________________________________________________________________
const TCollection*
AliMUONTOTCAStoreIterator::GetCollection() const
{
  /// The top level collection we're iterating upon, i.e. a TObjArray
  return fkData;
}

//_____________________________________________________________________________
TObject*
AliMUONTOTCAStoreIterator::Next()
{
  /// Find and return next element in the store
  
  if ( fCurrentTCA && fCurrentTCAIndex < fCurrentTCA->GetLast() ) 
  {
    ++fCurrentTCAIndex;
  }
  else
  {
    fCurrentTCAIndex = 0;
    fCurrentTCA = 0;
    
    while ( ( !fCurrentTCA || fCurrentTCA->GetLast()==-1 ) && 
            fCurrentChamberId < fLastChamberId ) 
    {
      ++fCurrentChamberId;
      fCurrentTCA = static_cast<TClonesArray*>(fkData->At(fCurrentChamberId));
    }
  }
  
  if ( fCurrentTCA ) 
  {
    // get the pointer to be returned
    return fCurrentTCA->At(fCurrentTCAIndex);
  }
  
  return 0x0;
}

//_____________________________________________________________________________
void
AliMUONTOTCAStoreIterator::Reset()
{
  /// Reset the iterator
  fCurrentTCAIndex = -1;
  fCurrentChamberId = fFirstChamberId-1;
  fCurrentTCA = 0x0;
}
