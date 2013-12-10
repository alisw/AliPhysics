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
    aodTrack->ResetStatus(AliVTrack::kEMCALmatch); //MV: necessary?
    if(aodTrack->Pt()<fMinPtCut) 
      continue;
    if(aodTrack->GetTrackPtOnEMCal()>0) 
      continue;

    Double_t phi = aodTrack->Phi()*TMath::RadToDeg();
    if (TMath::Abs(aodTrack->Eta())>0.9 || phi <= 10 || phi >= 250) 
      continue;

    Double_t xyz[3], pxpypz[3], cv[21];
    aodTrack->PxPyPz(pxpypz);  
    aodTrack->XvYvZv(xyz);
    aodTrack->GetCovarianceXYZPxPyPz(cv);  
    AliExternalTrackParam *trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,aodTrack->Charge());
    //    AliExternalTrackParam *trackParam =  const_cast<AliExternalTrackParam*>(eTrack->GetInnerParam()); MV: note, not taking InnerParam in AOD, not available
    if(!trackParam) 
      continue;

    // Extrapolate the track to EMCal surface
    AliExternalTrackParam emcalParam(*trackParam);
    Float_t etaout=-999, phiout=-999, ptout=-999;
    Bool_t ret = fRecoUtils->ExtrapolateTrackToEMCalSurface(&emcalParam, 
                                                            fDist, 
                                                            aodTrack->M(), 
                                                            fRecoUtils->GetStepSurface(), 
                                                            etaout, 
                                                            phiout,
							    ptout);
    if (!ret)
      continue;
    if (TMath::Abs(etaout)>0.75 || (phiout<70*TMath::DegToRad()) || (phiout>190*TMath::DegToRad()))
      continue;
    aodTrack->SetTrackPhiEtaPtOnEMCal(phiout, etaout, ptout);
    aodTrack->SetStatus(AliVTrack::kEMCALmatch);

    delete trackParam;
  }
}
