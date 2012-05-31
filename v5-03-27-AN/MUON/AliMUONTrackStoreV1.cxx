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
/// \class AliMUONTrackStoreV1
///
/// Implementation of AliMUONTrackStoreV1, which should be backward
/// compatible, i.e. able to read old TreeT files
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONTrackStoreV1.h"

#include <TClonesArray.h>
#include <TTree.h>
#include "AliLog.h"
#include "AliMUONTrack.h"
#include "AliMUONTreeManager.h"

/// \cond CLASSIMP
ClassImp(AliMUONTrackStoreV1)
/// \endcond

//_____________________________________________________________________________
AliMUONTrackStoreV1::AliMUONTrackStoreV1() : AliMUONVTrackStore(),
 fTracks(0x0)
{
   /// Ctor
  CreateTracks();
}

//_____________________________________________________________________________
AliMUONTrackStoreV1::AliMUONTrackStoreV1(TRootIOCtor* /*dummy*/) : AliMUONVTrackStore(),
fTracks(0x0)
{
  /// Ctor
}

//_____________________________________________________________________________
AliMUONTrackStoreV1::~AliMUONTrackStoreV1()
{
  /// dtor
  delete fTracks;
}

//_____________________________________________________________________________
AliMUONTrack* 
AliMUONTrackStoreV1::Add(const AliMUONTrack& track)
{
  /// Add a track
  
  if (!fTracks) CreateTracks();

  return new((*fTracks)[fTracks->GetLast()+1]) AliMUONTrack(track);
}

//_____________________________________________________________________________
AliMUONTrack*
AliMUONTrackStoreV1::Remove(AliMUONTrack& track)
{
  /// Remove a track from the store
  AliMUONTrack* t = static_cast<AliMUONTrack*>(fTracks->Remove(&track));
  if (t) fTracks->Compress();
  return t;
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackStoreV1::Connect(TTree& tree, Bool_t alone) const
{
  /// Connect this store to the tree
  AliMUONTreeManager tman;
  
  Bool_t ok;
  
  if ( tree.GetBranch("MUONTrack") )
  {
    if ( alone ) tman.UpdateBranchStatuses(tree,"MUONTrack");
    ok = tman.SetAddress(tree,"MUONTrack",TracksPtr());
  }
  else
  {
    ok = tman.MakeBranch(tree,ClassName(),"TClonesArray","MUONTrack",
                         TracksPtr());
  }
  return ok;
}

//_____________________________________________________________________________
TIterator*
AliMUONTrackStoreV1::CreateIterator() const
{
  /// Create an iterator to loop over tracks
  if ( fTracks ) return fTracks->MakeIterator();
  return 0x0;
}

//_____________________________________________________________________________
void 
AliMUONTrackStoreV1::Clear(Option_t*)
{
  /// Reset
  if (fTracks) fTracks->Clear("C");
}

//_____________________________________________________________________________
void
AliMUONTrackStoreV1::CreateTracks()
{
  /// Allocate track container
  if (fTracks) 
  {
    AliError("Cannot allocate again fTracks as it is there already !");
  }
  else
  {
    fTracks = new TClonesArray("AliMUONTrack",10);
  }
}

//_____________________________________________________________________________
Int_t
AliMUONTrackStoreV1::GetSize() const
{
  /// Return the number of tracks we hold
  if ( fTracks ) return fTracks->GetLast()+1;
  return 0;
}
