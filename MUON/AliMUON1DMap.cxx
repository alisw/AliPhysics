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

#include "AliMUON1DMap.h"

#include "AliLog.h"
#include "AliMpExMap.h"
#include "AliMUON1DMapIterator.h"

///
/// \class AliMUON1DMap
/// This class is simply a wrapper to an AliMpExMap, offering in addition a
/// control over the replacement policy when you add
/// something to it.
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUON1DMap)
/// \endcond

//_____________________________________________________________________________
AliMUON1DMap::AliMUON1DMap(Int_t theSize)
: AliMUONV1DStore(),
  fMap(new AliMpExMap(true))
{
/// Default ctor

  if ( theSize ) 
  {
    fMap->SetSize(theSize);
  }
    fMap->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUON1DMap::AliMUON1DMap(const AliMUON1DMap& other)
: AliMUONV1DStore(),
  fMap(0x0)
{
/// Copy constructor

  other.CopyTo(*this);
}

//_____________________________________________________________________________
AliMUON1DMap&
AliMUON1DMap::operator=(const AliMUON1DMap& other)
{
/// Assignment operator

  other.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUON1DMap::~AliMUON1DMap()
{
/// dtor, we're the owner of our internal map.

  delete fMap;
}

//_____________________________________________________________________________
void
AliMUON1DMap::CopyTo(AliMUON1DMap& dest) const
{
/// Make a deep copy

  delete dest.fMap;
  dest.fMap = fMap;
}

//_____________________________________________________________________________
TObject* 
AliMUON1DMap::Get(Int_t i) const
{
/// Get the object located at index i, if it exists, and if i is correct.

  return fMap->GetValue(i);
}

//_____________________________________________________________________________
AliMUONVDataIterator*
AliMUON1DMap::Iterator() const
{
  // Create and return an iterator on this map
  // Returned iterator must be deleted by user.
  if ( fMap ) 
  {
    return new AliMUON1DMapIterator(*fMap);
  }
  return 0x0;
}


//_____________________________________________________________________________
Bool_t 
AliMUON1DMap::Set(Int_t i, TObject* object, Bool_t replace)
{
/// Set the object located at i
/// If replace=kFALSE and there's already an object at location i,
/// this method fails and returns kFALSE, otherwise it returns kTRUE

  TObject* o = Get(i);
  if ( o && !replace )
  {
    AliError(Form("Object %p is already there for i=%d",o,i));
    return kFALSE;
  }
  if ( replace ) 
  {
    delete o;
  }
  fMap->Add(i,object);
  return kTRUE;
}



