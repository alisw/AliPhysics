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
#include "AliMpExMap.h"

/// 
/// Basic implementation of AliMUONV2DStore container using
/// AliMpExMap internally.
/// What we store is a "double" map : an AliMpExMap of AliMpExMaps
///

ClassImp(AliMUON2DMap)

//_____________________________________________________________________________
AliMUON2DMap::AliMUON2DMap() : AliMUONV2DStore(), fMap(new AliMpExMap(true))
{
  //
  // ctor.
  //
}

//_____________________________________________________________________________
AliMUON2DMap::AliMUON2DMap(const AliMUON2DMap& other)
: AliMUONV2DStore(),
fMap(0x0)
{
  other.CopyTo(*this);
}

//_____________________________________________________________________________
AliMUON2DMap&
AliMUON2DMap::operator=(const AliMUON2DMap& other)
{
  other.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUON2DMap::~AliMUON2DMap()
{
  //
  // dtor. we delete the map, which will delete the objects, as we're owner.
  //
  delete fMap;
}

//_____________________________________________________________________________
void
AliMUON2DMap::CopyTo(AliMUON2DMap&) const
{
  // 
  // Copy this into dest.
  //
  AliFatal("Implement me if needed");
}

//_____________________________________________________________________________
TObject* 
AliMUON2DMap::Get(Int_t i, Int_t j) const
{
  //
  // Return the value at position (i,j).
  //
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
void
AliMUON2DMap::Print(Option_t*) const
{
  //
  // Not implemented (yet?)
  //
}

//_____________________________________________________________________________
Bool_t 
AliMUON2DMap::Set(Int_t i, Int_t j, TObject* object, Bool_t replace)
{
  //
  // Set the object at position (i,j).
  // If replace==kTRUE, we don't care if there's an object there already,
  // otherwise we might refuse to set if the (i,j) location is already
  // filled (in which case we return kFALSE).
  
  TObject* o = fMap->GetValue(i);
  if ( !o )
  {
    AliMpExMap* m = new AliMpExMap(true);
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





