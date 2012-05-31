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
/// \class AliMUONVTrackStore
///
/// Base class of a track store
///
/// Note that objects stored are of concrete class AliMUONTrack, which 
/// might evolve to a virtual AliMUONVTrack for instance in the future...
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONVTrackStore.h"
#include "AliMUONTrack.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONVTrackStore)
/// \endcond

//_____________________________________________________________________________
AliMUONVTrackStore::AliMUONVTrackStore()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVTrackStore::~AliMUONVTrackStore()
{
  /// dtor
}

//_____________________________________________________________________________
Bool_t
AliMUONVTrackStore::Add(TObject* object)
{
  /// Add object, if of type AliMUONTrack
  if (object)
  {
    AliMUONTrack* t = dynamic_cast<AliMUONTrack*>(object);
    if (t)
    {
      Add(*t);
      return kTRUE;
    }
    else
    {
      AliError(Form("object not of expected AliMUONTrack type but %s",object->ClassName()));
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________
AliMUONVTrackStore*
AliMUONVTrackStore::Create(TTree& tree)
{
  /// Create a VTrackStore from the tree (if possible)
  return static_cast<AliMUONVTrackStore*>(AliMUONVStore::Create(tree,"MUONTrack"));
}
