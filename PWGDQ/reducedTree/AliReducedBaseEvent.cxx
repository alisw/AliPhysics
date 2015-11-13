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

ClassImp(AliReducedBaseEvent)

TClonesArray* AliReducedBaseEvent::fgTracks = 0;

//____________________________________________________________________________
AliReducedBaseEvent::AliReducedBaseEvent() :
  TObject(),
  fEventTag(0),
  fRunNo(0),
  fVtx(),
  fCentrality(),
  fCentQuality(0),
  fNtracks(),
  fTracks(0x0)
{
  //
  // Constructor
  //
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.;}
  for(Int_t i=0; i<7; ++i) fCentrality[i]=-1.;
  fNtracks[0]=0; fNtracks[1]=0;
}


//____________________________________________________________________________
AliReducedBaseEvent::AliReducedBaseEvent(const Char_t* /*name*/) :
  TObject(),
  fEventTag(0),
  fRunNo(0),
  fVtx(),
  fCentrality(),
  fCentQuality(0),
  fNtracks(),
  fTracks(0x0)
{
  //
  // Constructor
  //
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.;}
  for(Int_t i=0; i<4; ++i) fCentrality[i]=-1.;
  fNtracks[0]=0; fNtracks[1]=0;
  
  //if(!fgTracks) fgTracks = new TClonesArray("AliReducedBaseTrack", 100000);
  //fTracks = fgTracks;
}


//____________________________________________________________________________
AliReducedBaseEvent::~AliReducedBaseEvent()
{
  //
  // De-Constructor
  //
  //ClearEvent();
}

//_____________________________________________________________________________
void AliReducedBaseEvent::ClearEvent() {
  //
  // clear the event
  //
  if(fTracks) fTracks->Clear("C");
  fEventTag = 0;
  fRunNo = 0;
  fCentQuality = 0;
  for(Int_t i=0;i<7;++i) fCentrality[i] = -1.0;
  fNtracks[0] = 0; fNtracks[1] = 0;
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.;}
}
