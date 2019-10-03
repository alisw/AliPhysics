// $Id$
//
// Class to make PicoTracks from jet 4-vectors
//
// Author: M. Verweij

#include <TClonesArray.h>
#include <TRandom3.h>
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliEmcalPicoTrackFromJetMaker.h"

ClassImp(AliEmcalPicoTrackFromJetMaker)

//________________________________________________________________________
AliEmcalPicoTrackFromJetMaker::AliEmcalPicoTrackFromJetMaker() : 
  AliAnalysisTaskSE("AliEmcalPicoTrackFromJetMaker"),
  fTracksOutName("PicoTracksFromJets"),
  fJetsInName("tracks"),
  fJetsIn(0),
  fTracksOut(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalPicoTrackFromJetMaker::AliEmcalPicoTrackFromJetMaker(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksOutName("PicoTracksFromJets"),
  fJetsInName("tracks"),
  fJetsIn(0),
  fTracksOut(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalPicoTrackFromJetMaker::~AliEmcalPicoTrackFromJetMaker()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalPicoTrackFromJetMaker::UserCreateOutputObjects()
{
  // Create my user objects.

  fTracksOut = new TClonesArray("AliPicoTrack");
  fTracksOut->SetName(fTracksOutName);
}

//________________________________________________________________________
void AliEmcalPicoTrackFromJetMaker::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }

  // retrieve tracks from input.
  if (!fJetsIn) { 
    fJetsIn = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsInName));
    if (!fJetsIn) {
      AliError(Form("Could not retrieve jets %s!", fJetsInName.Data())); 
      return;
    }
    if (!fJetsIn->GetClass()->GetBaseClass("AliVParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fJetsInName.Data())); 
      return;
    }
  }

  // add tracks to event if not yet there
  fTracksOut->Delete();
  if (!(InputEvent()->FindListObject(fTracksOutName))) {
    InputEvent()->AddObject(fTracksOut);
  }

  // loop over tracks
  const Int_t Njets = fJetsIn->GetEntriesFast();
  for (Int_t iJets = 0, nacc = 0; iJets < Njets; ++iJets) {

    AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJetsIn->At(iJets));
    if (!jet)
      continue;

    Bool_t isEmc = kFALSE;
    if (TMath::Abs(jet->Eta()) < 0.75 && 
	jet->Phi() > 70 * TMath::DegToRad() &&jet->Phi() < 190 * TMath::DegToRad())
      isEmc = kTRUE;

    AliPicoTrack *picotrack = new ((*fTracksOut)[nacc]) AliPicoTrack(jet->Pt(), 
								     jet->Eta(), 
								     jet->Phi(), 
								     1, 
								     1,
								     0,
								     jet->Eta(), 
								     jet->Phi(), 
								     jet->Pt(), 
								     isEmc,
								     jet->M());
    picotrack->SetTrackType(0);
    ++nacc;
  }
}
