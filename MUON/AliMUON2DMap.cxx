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

#include "AliMUON2DMap.h"

#include "AliLog.h"
#include "AliMUON2DMapIterator.h"
#include "AliMUON2DMapIteratorByI.h"
#include "AliMpExMap.h"

//-----------------------------------------------------------------------------
/// \class AliMUON2DMap
/// Basic implementation of AliMUONVStore container using
/// AliMpExMap internally.
/// What we store is a "double" map : an AliMpExMap of AliMpExMaps
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUON2DMap)
/// \endcond

namespace
{
  //___________________________________________________________________________
  TObject* GetValue(TExMapIter& iter, Int_t& theKey) 
  {
    /// return the next value corresponding to theKey in iterator iter
    theKey = -1;
    Long_t key, value;
    Bool_t ok = iter.Next(key,value);
    if (!ok) return 0x0;
    theKey = (Int_t)(key & 0xFFFF);
    return reinterpret_cast<TObject*>(value);
  }
}

//_____________________________________________________________________________
AliMUON2DMap::AliMUON2DMap(Bool_t optimizeForDEManu) 
: AliMUONVStore(), 
  fMap(0x0),
  fOptimizeForDEManu(optimizeForDEManu)
{
  /// Default constructor.
    Clear();
}

//_____________________________________________________________________________
AliMUON2DMap::AliMUON2DMap(const AliMUON2DMap& other)
: AliMUONVStore(),
fMap(0x0),
fOptimizeForDEManu(kFALSE)
{
 /// Copy constructor.

 other.CopyTo(*this);
}

//_____________________________________________________________________________
AliMUON2DMap&
AliMUON2DMap::operator=(const AliMUON2DMap& other)
{
/// Assignment operator

  other.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUON2DMap::~AliMUON2DMap()
{
/// Destructor. 
/// We delete the map, which will delete the objects, as we're owner.

  delete fMap;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUON2DMap::Create() const
{
  /// Create a void copy of *this. 
  return new AliMUON2DMap(fOptimizeForDEManu);
}

//_____________________________________________________________________________
void
AliMUON2DMap::CopyTo(AliMUON2DMap& dest) const
{
  /// Copy this into dest.

  delete dest.fMap;
  dest.fMap = new AliMpExMap(*fMap);
  dest.fOptimizeForDEManu = fOptimizeForDEManu;
}

//_____________________________________________________________________________
Bool_t
AliMUON2DMap::Add(TObject* object)
{
  /// Add object, using the decoding of uniqueID into two ints as the key
  UInt_t uniqueID = object->GetUniqueID();
  Int_t j = ( uniqueID & 0xFFFF0000 ) >> 16;
  Int_t i = ( uniqueID & 0xFFFF);
  return Set(i,j,object,kFALSE);
}

//_____________________________________________________________________________
TObject* 
AliMUON2DMap::FindObject(Int_t i, Int_t j) const
{
  /// Return the value at position (i,j).

  AliMpExMap* m = static_cast<AliMpExMap*>(fMap->GetValue(i));
  if (m) return m->GetValue(j);
  return 0x0;
}

//_____________________________________________________________________________
TIterator*
AliMUON2DMap::CreateIterator() const
{
  // Create and return an iterator on this map
  // Returned iterator must be deleted by user.
  if ( fMap ) 
  {
    return new AliMUON2DMapIterator(fMap);
  }
  return 0x0;
}

//_____________________________________________________________________________
TIterator*
AliMUON2DMap::CreateIterator(Int_t firstI, Int_t lastI) const
{
  // Create and return an iterator on this map
  // Returned iterator must be deleted by user.
  if ( fMap ) 
  {
    return new AliMUON2DMapIteratorByI(*fMap,firstI,lastI);
  }
  return 0x0;
}

//_____________________________________________________________________________
void 
AliMUON2DMap::Clear(Option_t*)
{
  /// Reset
  delete fMap;

  fMap = new AliMpExMap(kTRUE);

  if ( fOptimizeForDEManu )
  {
    fMap->SetSize(228); // hard-coded constant in order not to depend on mapping
    // if this number ever change, it will not break the code, simply the
    // automatic resizing will give a warning...
  }
}  

//_____________________________________________________________________________
Int_t 
AliMUON2DMap::GetSize() const
{
  /// Return the number of objects we hold
  TExMapIter iter(fMap->GetIterator());
  Int_t i;
  Int_t theSize(0);
  
  while ( GetValue(iter,i) ) 
  {
    theSize += GetSize(i);
  }
  return theSize;
}

//_____________________________________________________________________________
Int_t 
AliMUON2DMap::GetSize(Int_t i) const
{
  /// Return the number of objects we hold
  AliMpExMap* m = static_cast<AliMpExMap*>(fMap->GetValue(i));
  if (m)
  {
    return m->GetSize();
  }
  return 0;
}

//_____________________________________________________________________________
Bool_t 
AliMUON2DMap::Set(Int_t i, Int_t j, TObject* object, Bool_t replace)
{
/// Set the object at position (i,j).
/// If replace==kTRUE, we don't care if there's an object there already,
/// otherwise we might refuse to set if the (i,j) location is already
/// filled (in which case we return kFALSE).
  
  TObject* o = fMap->GetValue(i);
  if ( !o )
  {
    AliMpExMap* m = new AliMpExMap(true);
    if ( fOptimizeForDEManu ) 
    {
      m->SetSize(451); // same remark as for the SetSize in ctor...
    }
    fMap->Add(i,m);
    o = fMap->GetValue(i);
  }
  AliMpExMap* m = static_cast<AliMpExMap*>(o);
 
  o = m->GetValue(j);
  
  if ( !o )
  {
    m->Add(j,object);
  }
  else 
  {
    if ( replace ) 
    {
      delete o;
      m->Add(j,object);
    }
    else
    {
      AliError(Form("Object %p is already there for (i,j)=(%d,%d)",o,i,j));
      return kFALSE;
    }
  }

  return kTRUE;
}





