#include "AliToyMCEvent.h"

ClassImp(AliToyMCEvent);
Int_t AliToyMCEvent::evCounter = 0;

AliToyMCEvent::AliToyMCEvent()
  :TObject()
  ,fEventNumber(0)
  ,fT0(-1.)
  ,fX(-1000.)
  ,fY(-1000.)
  ,fZ(-1000.)
  ,fTracks("AliToyMCTrack")
{
  fEventNumber = evCounter;
  evCounter++;
}

//____________________________________________________
AliToyMCEvent::AliToyMCEvent(const AliToyMCEvent &event)
  : TObject(event)
  ,fEventNumber(event.fEventNumber)
  ,fT0(event.fT0)
  ,fX(event.fX)
  ,fY(event.fY)
  ,fZ(event.fZ)
  ,fTracks(event.fTracks)
{
  //
}

//_____________________________________________________
AliToyMCEvent& AliToyMCEvent::operator = (const AliToyMCEvent &event)
{
  //assignment operator
  if (&event == this) return *this;
  new (this) AliToyMCEvent(event);

  return *this;
}
//_____________________________________________________
AliToyMCTrack* AliToyMCEvent::AddTrack(const AliToyMCTrack &track)
{
  return new(fTracks[fTracks.GetEntriesFast()]) AliToyMCTrack(track);
}

