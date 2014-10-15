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
/// \class AliMUONVHitStore
///
/// Base class of a MUON hit store
///
/// \author Laurent Aphecetche, Subatech
///
//-----------------------------------------------------------------------------

#include "AliMUONVHitStore.h"
#include "AliMUONHit.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONVHitStore)
/// \endcond

//_____________________________________________________________________________
AliMUONVHitStore::AliMUONVHitStore()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVHitStore::~AliMUONVHitStore()
{
  /// dtor
}

//_____________________________________________________________________________
Bool_t
AliMUONVHitStore::Add(TObject* object)
{
  /// Add an object, if of the right type
  if (object)
  {
    AliMUONHit* hit = dynamic_cast<AliMUONHit*>(object);
    if (hit)
    {
      Add(*hit);
      return kTRUE;
    }
    else
    {
      AliError(Form("object not of expected AliMUONHit type but %s",object->ClassName()));
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________
AliMUONVHitStore*
AliMUONVHitStore::Create(TTree& tree)
{
  /// Create a VHitStore from the tree (if possible)
  return static_cast<AliMUONVHitStore*>(AliMUONVStore::Create(tree,"Hit"));
}
