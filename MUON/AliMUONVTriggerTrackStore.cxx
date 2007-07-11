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
/// \class AliMUONVTriggerTrackStore
///
/// Base class of a trigger track store
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONTriggerTrack.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONVTriggerTrackStore)
/// \endcond

//_____________________________________________________________________________
AliMUONVTriggerTrackStore::AliMUONVTriggerTrackStore()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVTriggerTrackStore::~AliMUONVTriggerTrackStore()
{
  /// dtor
}

//_____________________________________________________________________________
Bool_t
AliMUONVTriggerTrackStore::Add(TObject* object)
{
  /// Add an object, if it is of type AliMUONTriggerTrack
  if (object)
  {
    AliMUONTriggerTrack* tt = dynamic_cast<AliMUONTriggerTrack*>(object);
    if (tt)
    {
      Add(*tt);
      return kTRUE;
    }
    else
    {
      AliError(Form("object is not of expected AliMUONTriggerTrack type but %s",
                    object->ClassName()));
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________
AliMUONVTriggerTrackStore*
AliMUONVTriggerTrackStore::Create(TTree& tree)
{
  /// Create a VTriggerTrackStore from the tree (if possible)
  return static_cast<AliMUONVTriggerTrackStore*>(AliMUONVStore::Create(tree,"MUONTriggerTrack"));
}
