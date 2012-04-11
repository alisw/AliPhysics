// $Id$
//
// Task to propagate tracks to EMCAL surface.
//
//

#include <TClonesArray.h>
#include "AliAnalysisManager.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDEvent.h"
#include "AliEmcalTrackPropagatorTask.h"

ClassImp(AliEmcalTrackPropagatorTask)

//________________________________________________________________________
AliEmcalTrackPropagatorTask::AliEmcalTrackPropagatorTask() : 
  AliAnalysisTaskSE("AliEmcalTrackPropagatorTask"),
  fRecoUtils(0),
  fTracksName(),
  fDist(430),
  fMinPtCut(0.35),
  fEsdEv(0),
  fTracks(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalTrackPropagatorTask::AliEmcalTrackPropagatorTask(const char *name) : 
  AliAnalysisTaskSE("AliEmcalTrackPropagatorTask"),
  fRecoUtils(0),
  fTracksName("TpcSpdVertexConstrainedTracks"),
  fDist(430),
  fMinPtCut(0.35),
  fEsdEv(0),
  fTracks(0)
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

  delete fRecoUtils;
}

//________________________________________________________________________
void AliEmcalTrackPropagatorTask::UserCreateOutputObjects()
{
  // Create histograms.

  if (!fRecoUtils) {
    fRecoUtils = new AliEMCALRecoUtils;
    fRecoUtils->SetStep(25);
    AliInfo("No reco utils given, creating default utils");
  }
}

//________________________________________________________________________
void AliEmcalTrackPropagatorTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  fEsdEv = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fEsdEv) {
    AliError("Task works only on ESD events, returning");
    return;
  }

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }

  // get tracks from event if not yet there
  if (fTracksName == "Tracks")
    am->LoadBranch("Tracks");
  fTracks = dynamic_cast<TClonesArray*>((InputEvent()->FindListObject(fTracksName)));
  if (!fTracks) {
    AliError(Form("Could not get tracks %s, returning", fTracksName.Data()));
    return;
  }

  // Loop over all tracks
  const Int_t ntr = fTracks->GetEntries();
  for (Int_t i=0; i<ntr; ++i) {
    AliESDtrack *eTrack = static_cast<AliESDtrack*>(fTracks->At(i));
    if (!eTrack)
      continue;
    eTrack->ResetStatus(AliVTrack::kEMCALmatch);
    if(eTrack->Pt()<fMinPtCut) 
      continue;
    Double_t phi = eTrack->Phi()*TMath::RadToDeg();
    if (TMath::Abs(eTrack->Eta())>0.8 || phi <= 20 || phi >= 240) 
      continue;
    if (TMath::Abs(eTrack->Eta())>0.8 || phi <= 20 || phi >= 240) 
      continue;
    AliExternalTrackParam *trackParam =  const_cast<AliExternalTrackParam*>(eTrack->GetInnerParam());
    if(!trackParam) 
      continue;

    // Extrapolate the track to EMCal surface
    AliExternalTrackParam emcalParam(*trackParam);
    Float_t etaout=-999, phiout=-999;
    Bool_t ret = fRecoUtils->ExtrapolateTrackToEMCalSurface(&emcalParam, 
                                                            fDist, 
                                                            fRecoUtils->GetMass(), 
                                                            fRecoUtils->GetStepSurface(), 
                                                            etaout, 
                                                            phiout);
    if (!ret)
      continue;
    if (TMath::Abs(etaout)>0.75 || (phiout<70*TMath::DegToRad()) || (phiout>190*TMath::DegToRad()))
      continue;
    eTrack->SetOuterParam(&emcalParam,AliExternalTrackParam::kMultSec);
    eTrack->SetStatus(AliVTrack::kEMCALmatch);
  }
}
