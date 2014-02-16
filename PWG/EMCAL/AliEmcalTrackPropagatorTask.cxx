// $Id: AliEmcalTrackPropagatorTask.cxx | Mon Dec 9 12:59:28 2013 +0100 | Constantin Loizides  $
//
// Task to propagate tracks to EMCAL surface.
//
// Author: C.Loizides

#include "AliEmcalTrackPropagatorTask.h"
#include <TClonesArray.h>
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDEvent.h"

ClassImp(AliEmcalTrackPropagatorTask)

//________________________________________________________________________
AliEmcalTrackPropagatorTask::AliEmcalTrackPropagatorTask() : 
  AliAnalysisTaskSE("AliEmcalTrackPropagatorTask"),
  fTracksInName(),
  fTracksOutName(),
  fDist(440),
  fOnlyIfNotSet(kTRUE),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalTrackPropagatorTask::AliEmcalTrackPropagatorTask(const char *name) : 
  AliAnalysisTaskSE("AliEmcalTrackPropagatorTask"),
  fTracksInName(),
  fTracksOutName(),
  fDist(440),
  fOnlyIfNotSet(kTRUE),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.

  if (!name)
    return;

  SetName(name);

  fBranchNames = "ESD:AliESDHeader.,Tracks";
}

//________________________________________________________________________
AliEmcalTrackPropagatorTask::~AliEmcalTrackPropagatorTask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalTrackPropagatorTask::UserCreateOutputObjects()
{
  // User create output objects.
}

//________________________________________________________________________
void AliEmcalTrackPropagatorTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return;
    }
  }

  if (fTracksInName.Length()==0) {
    if (esdev) {
      fTracksInName = "Tracks";
    } else {
      fTracksInName = "tracks";
    }
  }

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }
  if (fTracksInName == "Tracks")
    am->LoadBranch("Tracks");
  else if (fTracksInName == "tracks")
    am->LoadBranch("tracks");

  if (!fTracksIn) {
    fTracksIn = dynamic_cast<TClonesArray*>((InputEvent()->FindListObject(fTracksInName)));
    if (!fTracksIn) {
      AliError(Form("Could not get tracks %s, returning", fTracksInName.Data()));
      return;
    }
  }

  if ((fTracksOutName.Length()>0) && !fTracksOut) {
    if ((InputEvent()->FindListObject(fTracksOutName))) {
      AliError(Form("Could not add tracks %s to event, returning", fTracksOutName.Data()));
      return;
    }
    if (esdev)
      fTracksOut = new TClonesArray("AliESDtrack");
    else 
      fTracksOut = new TClonesArray("AliAODTrack");
    fTracksOut->SetName(fTracksOutName);
    InputEvent()->AddObject(fTracksOut);
  }

  const Int_t ntr = fTracksIn->GetEntries();
  for (Int_t i=0; i<ntr; ++i) {
    AliVTrack *inTrack = dynamic_cast<AliVTrack*>(fTracksIn->At(i));
    if (!inTrack)
      continue;
    if (inTrack->IsExtrapolatedToEMCAL() && fOnlyIfNotSet)
      continue;
    AliVTrack *outTrack = inTrack;
    if (fTracksOut) {
      if (esdev)
	outTrack = new ((*fTracksOut)[i]) AliESDtrack(*((AliESDtrack*)inTrack));
      else
	outTrack = new ((*fTracksOut)[i]) AliAODTrack(*((AliAODTrack*)inTrack));
    }
    AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(outTrack,fDist);
  }
}
