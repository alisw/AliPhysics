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
/// \class AliMUONTriggerTrackStoreV1
///
/// Implementation of AliMUONVTriggerTrackStore which should be
/// backward compatible, i.e. able to read TreeT produced before
/// the introduction of the AliMUONVStore concept
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONTriggerTrackStoreV1.h"

#include <TClonesArray.h>
#include "AliMUONTreeManager.h"
#include "AliMUONTriggerTrack.h"
#include <TTree.h>

/// \cond CLASSIMP
ClassImp(AliMUONTriggerTrackStoreV1)
/// \endcond

//_____________________________________________________________________________
AliMUONTriggerTrackStoreV1::AliMUONTriggerTrackStoreV1(TRootIOCtor* /*dummy*/) : AliMUONVTriggerTrackStore(),
fTracks(0x0)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONTriggerTrackStoreV1::AliMUONTriggerTrackStoreV1() : AliMUONVTriggerTrackStore(),
 fTracks(new TClonesArray("AliMUONTriggerTrack",10))
{
   /// ctor
   fTracks->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONTriggerTrackStoreV1::~AliMUONTriggerTrackStoreV1()
{
  /// Dtor
  delete fTracks;
}

//_____________________________________________________________________________
void 
AliMUONTriggerTrackStoreV1::Add(const AliMUONTriggerTrack& track)
{
  /// Add a new trigger track
  new((*fTracks)[fTracks->GetLast()+1]) AliMUONTriggerTrack(track);
}

//_____________________________________________________________________________
Bool_t
AliMUONTriggerTrackStoreV1::Connect(TTree& tree, Bool_t alone) const
{
  /// Connect this to the tree
  AliMUONTreeManager tman;
  Bool_t ok;
  
  if ( tree.GetBranch("MUONTriggerTrack") ) 
  {
    if ( alone ) tman.UpdateBranchStatuses(tree,"MUONTriggerTrack");
    ok = tman.SetAddress(tree,"MUONTriggerTrack",TracksPtr());
  }
  else
  {
    ok = tman.MakeBranch(tree,ClassName(),"TClonesArray","MUONTriggerTrack",
                         TracksPtr());
  }

  return kTRUE;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerTrackStoreV1::GetSize() const
{
  /// Return the number of trigger tracks we hold
  return fTracks->GetLast()+1;
}

//_____________________________________________________________________________
TIterator*
AliMUONTriggerTrackStoreV1::CreateIterator() const
{
  /// Return an iterator to loop over trigger tracks
  return fTracks->MakeIterator();
}

//_____________________________________________________________________________
void 
AliMUONTriggerTrackStoreV1::Clear(Option_t*)
{
  /// Reset
  fTracks->Clear("C");
}
