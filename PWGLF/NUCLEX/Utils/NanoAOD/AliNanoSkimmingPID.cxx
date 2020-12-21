#include "AliNanoSkimmingPID.h"

#include <TObject.h>

#include <AliVEvent.h>

AliNanoSkimmingPID::AliNanoSkimmingPID()
    : AliAnalysisCuts{}, fTrackFilter{} {

}

bool AliNanoSkimmingPID::IsSelected(TObject *obj) {
  fSelected = false;
  fFilterMask = 0u;
  AliVEvent *ev = dynamic_cast<AliVEvent *>(obj);
  if (!ev) {
    return fSelected;
  }
  for (int iTrack{0}; iTrack < ev->GetNumberOfTracks() && !fSelected; ++iTrack)
    fSelected = fTrackFilter.IsSelected(ev->GetTrack(iTrack));

  return fSelected;
}

bool AliNanoSkimmingPID::IsSelected(TList *) {
  Fatal("AliNanoSkimmingPID::IsSelected","Method not implemented for lists");
  return false;
}
