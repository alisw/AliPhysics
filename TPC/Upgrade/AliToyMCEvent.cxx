#include "AliToyMCEvent.h"

ClassImp(AliToyMCEvent);
Int_t AliToyMCEvent::fgEvCounter = 0;

AliToyMCEvent::AliToyMCEvent()
  :TObject()
  ,fEventNumber(0)
  ,fEventType(kPhysics)
  ,fT0(-1.)
  ,fX(-1000.)
  ,fY(-1000.)
  ,fZ(-1000.)
  ,fSCscale(-1.)
  ,fSCscaleChi2(0)
  ,fTracks("AliToyMCTrack")
{
  fEventNumber = fgEvCounter;
  fgEvCounter++;
}

//____________________________________________________
AliToyMCEvent::AliToyMCEvent(const AliToyMCEvent &event)
  : TObject(event)
  ,fEventNumber(event.fEventNumber)
  ,fEventType(event.fEventType)
  ,fT0(event.fT0)
  ,fX(event.fX)
  ,fY(event.fY)
  ,fZ(event.fZ)
  ,fSCscale(event.fSCscale)
  ,fSCscaleChi2(event.fSCscaleChi2)
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
//____________________________________________________
AliToyMCTrack* AliToyMCEvent::AddTrack(Double_t xyz[3],Double_t pxpypz[3],
                        Double_t cv[21],Short_t sign)
{
  return new(fTracks[fTracks.GetEntriesFast()]) AliToyMCTrack(xyz,pxpypz,cv,sign);
}

