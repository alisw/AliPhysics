#include "ToyMCEvent.h"

ClassImp(ToyMCEvent);
Int_t ToyMCEvent::evCounter = 0;

ToyMCEvent::ToyMCEvent()
  :TObject()
  ,fT0(-1.)
  ,fX(-1000.)
  ,fY(-1000.)
  ,fZ(-1000.)
  ,fTracks("ToyMCTrack")
{
  fEventNumber = evCounter;
  evCounter++;
}

//____________________________________________________
ToyMCEvent::ToyMCEvent(const ToyMCEvent &event)
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
ToyMCEvent& ToyMCEvent::operator = (const ToyMCEvent &event)
{
  //assignment operator
  if (&event == this) return *this;
  new (this) ToyMCEvent(event);

  return *this;
}
//_____________________________________________________
ToyMCTrack* ToyMCEvent::AddTrack(const ToyMCTrack &track)
{
  return new(fTracks[fTracks.GetEntriesFast()]) ToyMCTrack(track);
}

