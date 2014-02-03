// $Id$
//
// Class to make PicoTracks in AOD/ESD events.
//
// Author: S.Aiola, C.Loizides

#include <TClonesArray.h>
#include <TRandom3.h>
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalPicoTrackMaker.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVTrack.h"

ClassImp(AliEmcalPicoTrackMaker)

//________________________________________________________________________
AliEmcalPicoTrackMaker::AliEmcalPicoTrackMaker() : 
  AliAnalysisTaskSE("AliEmcalPicoTrackMaker"),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fMinTrackPt(0),
  fMaxTrackPt(1000),
  fMinTrackEta(-10),
  fMaxTrackEta(10),
  fMinTrackPhi(-10),
  fMaxTrackPhi(10),
  fTrackEfficiency(1),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.

  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
}

//________________________________________________________________________
AliEmcalPicoTrackMaker::AliEmcalPicoTrackMaker(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fMinTrackPt(0),
  fMaxTrackPt(1000),
  fMinTrackEta(-10),
  fMaxTrackEta(10),
  fMinTrackPhi(-10),
  fMaxTrackPhi(10),
  fTrackEfficiency(1),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.

  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
  fBranchNames = "ESD:AliESDHeader.,AliESDRun.,SPDVertex.,Tracks";
}

//________________________________________________________________________
AliEmcalPicoTrackMaker::~AliEmcalPicoTrackMaker()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalPicoTrackMaker::UserCreateOutputObjects()
{
  // Create my user objects.

  fTracksOut = new TClonesArray("AliPicoTrack");
  fTracksOut->SetName(fTracksOutName);
}

//________________________________________________________________________
void AliEmcalPicoTrackMaker::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }

  // retrieve tracks from input.
  if (!fTracksIn) { 
    fTracksIn = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksInName));
    if (!fTracksIn) {
      AliError(Form("Could not retrieve tracks %s!", fTracksInName.Data())); 
      return;
    }
    if (!fTracksIn->GetClass()->GetBaseClass("AliVParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fTracksInName.Data())); 
      return;
    }
  }

  // add tracks to event if not yet there
  fTracksOut->Delete();
  if (!(InputEvent()->FindListObject(fTracksOutName))) {
    InputEvent()->AddObject(fTracksOut);
  }

  // loop over tracks
  const Int_t Ntracks = fTracksIn->GetEntriesFast();
  for (Int_t iTracks = 0, nacc = 0; iTracks < Ntracks; ++iTracks) {

    AliVTrack *track = static_cast<AliVTrack*>(fTracksIn->At(iTracks));

    if (!track)
      continue;

    if (track->Pt() > fMaxTrackPt || track->Pt() < fMinTrackPt)
      continue;

    if (track->Eta() < fMinTrackEta || track->Eta() > fMaxTrackEta || 
	track->Phi() < fMinTrackPhi || track->Phi() > fMaxTrackPhi)
      continue;

    if (fTrackEfficiency < 1) {
      Double_t r = gRandom->Rndm();
      if (fTrackEfficiency < r) 
	continue;
    }

    Bool_t isEmc = kFALSE;
    if (TMath::Abs(track->GetTrackEtaOnEMCal()) < 0.75 && 
	track->GetTrackPhiOnEMCal() > 70 * TMath::DegToRad() &&
	track->GetTrackPhiOnEMCal() < 190 * TMath::DegToRad())
      isEmc = kTRUE;

    AliPicoTrack *picotrack = new ((*fTracksOut)[nacc]) AliPicoTrack(track->Pt(), 
								     track->Eta(), 
								     track->Phi(), 
								     track->Charge(), 
								     track->GetLabel(),
								     AliPicoTrack::GetTrackType(track),
								     track->GetTrackEtaOnEMCal(), 
								     track->GetTrackPhiOnEMCal(), 
								     track->GetTrackPtOnEMCal(), 
								     isEmc);
    picotrack->SetTrack(track);
    ++nacc;
  }
}
