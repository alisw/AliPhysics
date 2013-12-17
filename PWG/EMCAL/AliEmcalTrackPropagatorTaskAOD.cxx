// $Id$
//
// Task to propagate AOD tracks to EMCAL surface.
//
// Author: C.Loizides

#include <TClonesArray.h>
#include "AliAnalysisManager.h"
#include "AliEMCALRecoUtils.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliExternalTrackParam.h"
#include <AliMagF.h>
#include <AliTrackerBase.h>
#include <AliEMCALRecoUtils.h>

#include "AliEmcalTrackPropagatorTaskAOD.h"

ClassImp(AliEmcalTrackPropagatorTaskAOD)

//________________________________________________________________________
AliEmcalTrackPropagatorTaskAOD::AliEmcalTrackPropagatorTaskAOD() : 
  AliAnalysisTaskSE("AliEmcalTrackPropagatorTaskAOD"),
  fRecoUtils(0),
  fTracksName(),
  fDist(430),
  fMinPtCut(0.35),
  fAodEv(0),
  fTracks(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalTrackPropagatorTaskAOD::AliEmcalTrackPropagatorTaskAOD(const char *name) : 
  AliAnalysisTaskSE("AliEmcalTrackPropagatorTaskAOD"),
  fRecoUtils(0),
  fTracksName("TpcSpdVertexConstrainedTracks"),
  fDist(430),
  fMinPtCut(0.35),
  fAodEv(0),
  fTracks(0)
{
  // Constructor.

  if (!name)
    return;

  SetName(name);

  //  fBranchNames = "ESD:AliESDHeader.,Tracks";
}

//________________________________________________________________________
AliEmcalTrackPropagatorTaskAOD::~AliEmcalTrackPropagatorTaskAOD()
{
  // Destructor.

  delete fRecoUtils;
}

//________________________________________________________________________
void AliEmcalTrackPropagatorTaskAOD::UserCreateOutputObjects()
{
  // Create histograms.

  if (!fRecoUtils) {
    fRecoUtils = new AliEMCALRecoUtils;
    fRecoUtils->SetStep(20);
    AliInfo("No reco utils given, creating default utils");
  }
}

//________________________________________________________________________
void AliEmcalTrackPropagatorTaskAOD::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAodEv) {
    AliError("Task works only on AOD events, returning");
    return;
  }

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }

  // get tracks from event if not yet there
  if (fTracksName == "tracks")
    am->LoadBranch("tracks");

  fTracks = dynamic_cast<TClonesArray*>((InputEvent()->FindListObject(fTracksName)));
  if (!fTracks) {
    AliError(Form("Could not get tracks %s, returning", fTracksName.Data()));
    return;
  }

  // Loop over all tracks
  const Int_t ntr = fTracks->GetEntries();
  for (Int_t i=0; i<ntr; ++i) {
    AliAODTrack *aodTrack = static_cast<AliAODTrack*>(fTracks->At(i));
    if (!aodTrack)
      continue;
    if(aodTrack->Pt()<fMinPtCut || aodTrack->GetTrackPtOnEMCal()>0) 
      continue;

    AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(aodTrack,fDist);
  }
}
