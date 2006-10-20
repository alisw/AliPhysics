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
#include "AliMUONVDataIterator.h"
#include "AliMUON2DMapIterator.h"
#include "AliMpExMap.h"
#include "AliMpIntPair.h"
#include "AliMpManuList.h"
#include "AliMpDEManager.h"
#include "AliMUONConstants.h"

/// \class AliMUON2DMap
/// \brief Basic implementation of AliMUONV2DStore container using
/// AliMpExMap internally.
/// What we store is a "double" map : an AliMpExMap of AliMpExMaps
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUON2DMap)
/// \endcond

//_____________________________________________________________________________
AliMUON2DMap::AliMUON2DMap(Bool_t optimizeForDEManu) 
: AliMUONV2DStore(), 
  fMap(new AliMpExMap(true)), 
  fOptimizeForDEManu(optimizeForDEManu)
{
/// Default constructor.
    if ( fOptimizeForDEManu )
    {
      Int_t nDEs(0);
      for ( Int_t i = 0; i < AliMUONConstants::NTrackingCh(); ++i )
      {
        nDEs += AliMpDEManager::GetNofDEInChamber(i);
      }
      fMap->SetSize(nDEs);
    }
}

//_____________________________________________________________________________
AliMUON2DMap::AliMUON2DMap(const AliMUON2DMap& other)
: AliMUONV2DStore(),
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
AliMUONV2DStore*
AliMUON2DMap::CloneEmpty() const
{
  /// Create a void copy of *this. 
  return new AliMUON2DMap;
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
TObject* 
AliMUON2DMap::Get(Int_t i, Int_t j) const
{
/// Return the value at position (i,j).

  TObject* o = fMap->GetValue(i);
  if ( o )
  {
    AliMpExMap* m = dynamic_cast<AliMpExMap*>(o);
    if (!m) AliFatal(Form("fMap[%d] not of the expected type",i));
    return m->GetValue(j);
  }
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVDataIterator*
AliMUON2DMap::Iterator() const
{
  // Create and return an iterator on this map
  // Returned iterator must be deleted by user.
  if ( fMap ) 
  {
    return new AliMUON2DMapIterator(*fMap);
  }
  return 0x0;
}

//_____________________________________________________________________________
void
AliMUON2DMap::Print(Option_t*) const
{
/// Not implemented (yet?)
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
      Int_t n(AliMpManuList::NumberOfManus(i));
      if (!n)
      {
        AliError(Form("This does not look right : i = %d is supposed to "
                      "be a DetElemId with n = %d manus!",i,n));
      }
      else
      {
        m->SetSize(n);
      }
    }
    fMap->Add(i,m);
    o = fMap->GetValue(i);
  }
  AliMpExMap* m = dynamic_cast<AliMpExMap*>(o);
  if (!m) AliFatal(Form("fMap[%d] not of the expected type",i));
  o = m->GetValue(j);
  if ( !o || ( o && replace ) )
  {
    if ( IsOwner() ) 
    {
      delete o;
    }
    m->Add(j,object);
  }
  else if ( o && !replace )
  {
    AliError(Form("Object %p is already there for (i,j)=(%d,%d)",o,i,j));
    return kFALSE;
  }
  return kTRUE;
}





