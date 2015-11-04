#ifndef ALIFILTEREDEVENT_H
#include "AliFilteredEvent.h"
#endif
#include <AliFilteredTrack.h>


ClassImp(AliFilteredEvent)

TClonesArray* AliFilteredEvent::fgTracks = 0;

AliFilteredEvent::AliFilteredEvent() :
  TObject(),
  fCentrality(-1.0),
  fVertex(),
  fRunNr(0),
  fEventPlaneAngle(-10.0),
  fNTracks(0),
  fTracks(0x0)
{
  //Constructor
  if(!fgTracks) fgTracks = new TClonesArray("AliFilteredTrack", 100000);
  fTracks = fgTracks;
  
}

AliFilteredEvent::AliFilteredEvent(const char* name) :
  TObject(),
  fCentrality(-1.0),
  fVertex(),
  fRunNr(0),
  fEventPlaneAngle(-10.0),
  fNTracks(0),
  fTracks(0x0)
{
  //Constructor
  if(!fgTracks) fgTracks = new TClonesArray("AliFilteredTrack", 100000);
  fTracks = fgTracks;
  
}
//____________________________________________________________________________
AliFilteredEvent::~AliFilteredEvent()
{
  //
  // De-Constructor
  //
  ClearEvent();
}


void AliFilteredEvent::ClearEvent() {
  fCentrality = -1.0;
  fVertex[0] = -20.0;
  fVertex[1] = -20.0;
  fVertex[3] = -20.0;
  fRunNr = 0;
  fEventPlaneAngle = -10.0;
  fNTracks = 0;
  if(fTracks) fTracks->Clear("C");
}