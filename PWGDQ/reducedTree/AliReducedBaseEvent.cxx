/*
***********************************************************
  Implementation of AliReducedBaseEvent class.
  Contact: iarsene@cern.ch
  2015/04/15
  *********************************************************
*/

#ifndef ALIREDUCEDBASEEVENT_H
#include "AliReducedBaseEvent.h"
#endif

#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"

ClassImp(AliReducedBaseEvent)

TClonesArray* AliReducedBaseEvent::fgTracks = 0;
TClonesArray* AliReducedBaseEvent::fgCandidates = 0;

//____________________________________________________________________________
AliReducedBaseEvent::AliReducedBaseEvent() :
  TObject(),
  fEventTag(0),
  fRunNo(0),
  fVtx(),
  fNVtxContributors(0),
  fCentrality(),
  fCentQuality(0),
  fNtracks(),
  fNV0candidates(),
  fTracks(0x0),
  fCandidates(0x0)
{
  //
  // Constructor
  //
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.;}
  for(Int_t i=0; i<7; ++i) fCentrality[i]=-1.;
  fNtracks[0]=0; fNtracks[1]=0;
  fNV0candidates[0]=0; fNV0candidates[1]=0;
}


//____________________________________________________________________________
AliReducedBaseEvent::AliReducedBaseEvent(const Char_t* /*name*/, Int_t trackOption /*=kNoInit*/) :
  TObject(),
  fEventTag(0),
  fRunNo(0),
  fVtx(),
  fNVtxContributors(0),
  fCentrality(),
  fCentQuality(0),
  fNtracks(),
  fNV0candidates(),
  fTracks(0x0),
  fCandidates(0x0)
{
  //
  // Constructor
  //
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.;}
  for(Int_t i=0; i<4; ++i) fCentrality[i]=-1.;
  fNtracks[0]=0; fNtracks[1]=0;
  fNV0candidates[0]=0; fNV0candidates[1]=0;
  if(trackOption == kUseBaseTracks) {
    if(!fgTracks) fgTracks = new TClonesArray("AliReducedBaseTrack", 100000);
    fTracks = fgTracks;
  }
  if(trackOption == kUseReducedTracks) {
     if(!fgTracks) fgTracks = new TClonesArray("AliReducedTrackInfo", 100000);
     fTracks = fgTracks;
  }
  if(!fgCandidates) fgCandidates = new TClonesArray("AliReducedPairInfo", 100000);
  fCandidates = fgCandidates;
}


//____________________________________________________________________________
AliReducedBaseEvent::~AliReducedBaseEvent()
{
  //
  // De-Constructor
  //
}

//_____________________________________________________________________________
void AliReducedBaseEvent::ClearEvent() {
  //
  // clear the event
  //
  if(fTracks) fTracks->Clear("C");
  if(fCandidates) fCandidates->Clear("C");
  fEventTag = 0;
  fRunNo = 0;
  fCentQuality = 0;
  for(Int_t i=0;i<7;++i) fCentrality[i] = -1.0;
  fNtracks[0] = 0; fNtracks[1] = 0;
  fNV0candidates[0] = 0; fNV0candidates[1] = 0;
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.;}
  fNVtxContributors = 0;
}
