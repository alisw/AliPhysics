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
TClonesArray* AliReducedBaseEvent::fgTracks2 = 0;
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
  fTracks2(0x0),
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
AliReducedBaseEvent::AliReducedBaseEvent(const Char_t* /*name*/, Int_t trackOption /*=kNoInit*/, Int_t track2Option /*=kNoInit*/) :
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
  fTracks2(0x0),
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
  
  if(track2Option == kUseBaseTracks || track2Option == kNoInit) {
     if(!fgTracks2) fgTracks2 = new TClonesArray("AliReducedBaseTrack", 100000);
     fTracks2 = fgTracks2;
  }
  if(track2Option == kUseReducedTracks) {
     if(!fgTracks2) fgTracks2 = new TClonesArray("AliReducedTrackInfo", 100000);
     fTracks2 = fgTracks2;
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

//____________________________________________________________________________
void AliReducedBaseEvent::CopyEventHeader(const AliReducedBaseEvent* other) {
   //
   // assignment operator overloading
   // NOTE: the fTracks and fCandidates arrays are not copied
   //
   ClearEvent();
   fEventTag = other->fEventTag;
   fRunNo = other->fRunNo;
   for(Int_t i=0; i<3; ++i) fVtx[i] = other->fVtx[i];
   fNVtxContributors = other->fNVtxContributors;
   for(Int_t i=0; i<7; ++i) fCentrality[i] = other->fCentrality[i];
   fCentQuality = other->fCentQuality;
   fNtracks[0] = other->fNtracks[0]; 
   fNV0candidates[0] = other->fNV0candidates[0];
}

//_____________________________________________________________________________
void AliReducedBaseEvent::ClearEvent() {
  //
  // clear the event
  //
  if(fTracks) fTracks->Clear("C");
  if(fTracks2) fTracks2->Clear("C");
  if(fCandidates) fCandidates->Clear("C");
  fEventTag = 0;
  fRunNo = 0;
  fCentQuality = 0;
  for(Int_t i=0;i<7;++i) fCentrality[i] = -9999.0;
  fNtracks[0] = 0; fNtracks[1] = 0;
  fNV0candidates[0] = 0; fNV0candidates[1] = 0;
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-9999.;}
  fNVtxContributors = 0;
}
