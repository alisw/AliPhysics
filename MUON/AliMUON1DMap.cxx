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
#include "AliMpExMap.h"
#include "AliMpExMapIterator.h"

#include "AliLog.h"

//-----------------------------------------------------------------------------
/// \class AliMUON1DMap
/// This class is simply a wrapper to an AliMpExMap, offering in addition a
/// control over the replacement policy when you add
/// something to it.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUON1DMap)
/// \endcond

//_____________________________________________________________________________
AliMUON1DMap::AliMUON1DMap(TRootIOCtor*)
: AliMUONVStore(),
fMap(0x0)
{
  /// I/O ctor
  
}

//_____________________________________________________________________________
AliMUON1DMap::AliMUON1DMap(Int_t theSize)
: AliMUONVStore(),
  fMap(new AliMpExMap)
{
/// Default ctor

  if ( theSize > 0) 
  {
    fMap->SetSize(theSize);
  }
  fMap->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUON1DMap::AliMUON1DMap(const AliMUON1DMap& other)
: AliMUONVStore(),
  fMap(new AliMpExMap(*other.fMap))
{
/// Copy constructor
}

//_____________________________________________________________________________
AliMUON1DMap&
AliMUON1DMap::operator=(const AliMUON1DMap& other)
{
/// Assignment operator
  *fMap = *other.fMap;
  return *this;
}

//_____________________________________________________________________________
AliMUON1DMap::~AliMUON1DMap()
{
/// destructor
  delete fMap;
}

//_____________________________________________________________________________
Bool_t 
AliMUON1DMap::Add(TObject* object)
{
  /// Add an object to this, using uniqueID as the key
  if (!object) return kFALSE;
  return Set(object->GetUniqueID(),object);
}

//_____________________________________________________________________________
void 
AliMUON1DMap::Clear(Option_t*)
{
  /// Reset
  fMap->Clear();
}

//_____________________________________________________________________________
AliMUON1DMap* 
AliMUON1DMap::Create() const
{
  /// Create an empty clone of this
  return new AliMUON1DMap(fMap->GetSize());
}

//_____________________________________________________________________________
TObject* 
AliMUON1DMap::FindObject(UInt_t i) const
{
/// Get the object located at index i, if it exists, and if i is correct.
  return fMap->GetValue(i);
}

//_____________________________________________________________________________
TObject* 
AliMUON1DMap::FindObject(Int_t i, Int_t j) const
{
  /// Get the object located at index (i,j), if it exists, and if i,j is correct.
  
  UInt_t uid = ( ( ( j & 0xFFFF ) << 16 ) | ( i & 0xFFFF ) );
  
  return fMap->GetValue(uid);
}

//_____________________________________________________________________________
TIterator*
AliMUON1DMap::CreateIterator() const
{
  /// Create and return an iterator on this map
  /// Returned iterator must be deleted by user.
  return fMap->CreateIterator();
}

//_____________________________________________________________________________
Int_t
AliMUON1DMap::GetSize() const
{
  /// Return the number of objects we hold
  return fMap->GetSize();
}

//_____________________________________________________________________________
Bool_t 
AliMUON1DMap::Set(Int_t i, TObject* object)
{
/// Set the object located at i
/// If there's already an object at location i,
/// this method fails and returns kFALSE, otherwise it returns kTRUE
  
  TObject* o = FindObject(i);
  if ( o )
  {
    AliError(Form("Object %p is already there for i=%d",o,i));
    return kFALSE;
  }
  fMap->Add(i,object);
  return kTRUE;
}

