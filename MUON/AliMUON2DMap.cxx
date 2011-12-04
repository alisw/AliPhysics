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
#include "AliMpExMapIterator.h"

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

const Int_t AliMUON2DMap::fgkOptimalSizeForDEManu = 228;

//_____________________________________________________________________________
AliMUON2DMap::AliMUON2DMap(TRootIOCtor*)
: AliMUONVStore(), 
fMap(0x0),
fOptimizeForDEManu(kFALSE)
{
  /// Root I/O constructor.
}

//_____________________________________________________________________________
AliMUON2DMap::AliMUON2DMap(Bool_t optimizeForDEManu) 
: AliMUONVStore(), 
  fMap(new AliMpExMap),
  fOptimizeForDEManu(optimizeForDEManu)
{
  /// Default constructor.
  // hard-coded constant in order not to depend on mapping
  // if this number ever change, it will not break the code, simply the
  // automatic resizing will give a warning...
    
  if ( fOptimizeForDEManu ) fMap->SetSize(fgkOptimalSizeForDEManu); 
}

//_____________________________________________________________________________
AliMUON2DMap::AliMUON2DMap(const AliMUON2DMap& other)
: AliMUONVStore(),
  fMap(new AliMpExMap(*other.fMap)),
  fOptimizeForDEManu(other.fOptimizeForDEManu)
{
 /// Copy constructor.
}

//_____________________________________________________________________________
AliMUON2DMap&
AliMUON2DMap::operator=(const AliMUON2DMap& other)
{
/// Assignment operator
  if ( this != &other )
  {
    *fMap = *other.fMap;
    fOptimizeForDEManu = other.fOptimizeForDEManu;
  }
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
Bool_t
AliMUON2DMap::Add(TObject* object)
{
  /// Add object, using the decoding of uniqueID into two ints as the key
  if (!object) return kFALSE;
  UInt_t uniqueID = object->GetUniqueID();
  Int_t j = ( uniqueID & 0xFFFF0000 ) >> 16;
  Int_t i = ( uniqueID & 0xFFFF);
  return Set(i,j,object,kFALSE);
}

//_____________________________________________________________________________
TObject* 
AliMUON2DMap::FindObject(UInt_t uid) const
{
  /// Return the value at position uid
  
  Int_t j = ( uid & 0xFFFF0000 ) >> 16;
  Int_t i = ( uid & 0xFFFF);
  return FindObject(i,j);
}

//_____________________________________________________________________________
TObject* 
AliMUON2DMap::FindObject(Int_t i, Int_t j) const
{
  /// Return the value at position (i,j).
  AliMpExMap* m = static_cast<AliMpExMap*>(fMap->GetValue(i));
  return m ? m->GetValue(j) : 0x0;
}

//_____________________________________________________________________________
TIterator*
AliMUON2DMap::CreateIterator() const
{
  // Create and return an iterator on this map
  // Returned iterator must be deleted by user.
  return new AliMUON2DMapIterator(*fMap);
}

//_____________________________________________________________________________
TIterator*
AliMUON2DMap::CreateIterator(Int_t firstI, Int_t lastI) const
{
  // Create and return an iterator on this map
  // Returned iterator must be deleted by user.
  return new AliMUON2DMapIteratorByI(*fMap,firstI,lastI);
}

//_____________________________________________________________________________
void 
AliMUON2DMap::Clear(Option_t*)
{
  /// Clear memory
  fMap->Clear();
}  

//_____________________________________________________________________________
Int_t 
AliMUON2DMap::GetSize() const
{
  /// Return the number of objects we hold
  TIter next(fMap->CreateIterator());
  Int_t theSize(0);
  AliMpExMap* m;
  
  while ( ( m = static_cast<AliMpExMap*>(next()) ) )
  {
    TIter next2(m->CreateIterator());
    while ( next2() ) 
    {
      ++theSize;
    }
  }
  return theSize;
}

//_____________________________________________________________________________
Int_t 
AliMUON2DMap::GetSize(Int_t i) const
{
  /// Return the number of objects we hold
  AliMpExMap* m = static_cast<AliMpExMap*>(fMap->GetValue(i));
  return m ? m->GetSize() : 0;
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
    AliMpExMap* m = new AliMpExMap;
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

