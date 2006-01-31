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

ClassImp(AliMUON1DMap)

//_____________________________________________________________________________
AliMUON1DMap::AliMUON1DMap() : AliMUONV1DStore(), fMap(new AliMpExMap(true))
{
}

//_____________________________________________________________________________
AliMUON1DMap::~AliMUON1DMap()
{
  delete fMap;
}

//_____________________________________________________________________________
TObject*
AliMUON1DMap::Get(Int_t detElemId) const
{
  return fMap->GetValue(detElemId);
}

//_____________________________________________________________________________
Bool_t
AliMUON1DMap::IsOwner() const
{
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t
AliMUON1DMap::Set(Int_t detElemId, TObject* object, Bool_t replace)
{
  TObject* o = Get(detElemId);
  if ( o && !replace )
  {
    AliError(Form("Object %p is already there for detElemId %d",o,detElemId));
    return kFALSE;
  }
  if ( o && IsOwner() ) 
  {
    delete o;
  }
  fMap->Add(detElemId,object);
  return kTRUE;
}





