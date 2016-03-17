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
#include "AliAODMCParticle.h"
#include "AliNamedArrayI.h"

ClassImp(AliEmcalPicoTrackMaker)

//________________________________________________________________________
AliEmcalPicoTrackMaker::AliEmcalPicoTrackMaker() : 
  AliAnalysisTaskSE("AliEmcalPicoTrackMaker"),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fMCParticlesName("mcparticles"),
  fMinTrackPt(0),
  fMaxTrackPt(1000),
  fMinTrackEta(-10),
  fMaxTrackEta(10),
  fMinTrackPhi(-10),
  fMaxTrackPhi(10),
  fTrackEfficiency(1),
  fCopyMCFlag(kFALSE),
  fTracksIn(0),
  fTracksOut(0),
  fMCParticles(0),
  fMCParticlesMap(0),
  fInit(kFALSE)
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
  fMCParticlesName("mcparticles"),
  fMinTrackPt(0),
  fMaxTrackPt(1000),
  fMinTrackEta(-10),
  fMaxTrackEta(10),
  fMinTrackPhi(-10),
  fMaxTrackPhi(10),
  fTrackEfficiency(1),
  fCopyMCFlag(kFALSE),
  fTracksIn(0),
  fTracksOut(0),
  fMCParticles(0),
  fMCParticlesMap(0),
  fInit(kFALSE)
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
}

//________________________________________________________________________
void AliEmcalPicoTrackMaker::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  if (!fInit) {
    fTracksIn = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksInName));
    if (!fTracksIn) {
      AliError(Form("Could not retrieve tracks %s!", fTracksInName.Data())); 
      return;
    }
    if (!fTracksIn->GetClass()->GetBaseClass("AliVParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fTracksInName.Data())); 
      return;
    }
    
    fTracksOut = new TClonesArray("AliPicoTrack");
    fTracksOut->SetName(fTracksOutName);    

    // add tracks to event if not yet there
    if (InputEvent()->FindListObject(fTracksOutName)) {
      AliFatal(Form("Object %s already present in the event!",fTracksOutName.Data()));
    }
    else {
      InputEvent()->AddObject(fTracksOut);
    }

    if (fCopyMCFlag) {
      fMCParticles = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCParticlesName));
      if (!fMCParticles) {
	AliError(Form("Could not retrieve MC particles %s!", fMCParticlesName.Data())); 
      }
      if (!fMCParticles->GetClass()->GetBaseClass("AliVParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fMCParticlesName.Data()));
	fMCParticles = 0;
      }
      
      TString mapName(fMCParticlesName);
      mapName += "_Map";
      fMCParticlesMap = dynamic_cast<AliNamedArrayI*>(InputEvent()->FindListObject(mapName));
    }

    fInit = kTRUE;
  }

  fTracksOut->Delete();

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
    
    if (fCopyMCFlag && track->GetLabel() != 0) {
      AliVParticle *mcpart = GetMCParticle(TMath::Abs(track->GetLabel()));
      if (mcpart) {
	UInt_t mcFlag = mcpart->GetFlag();	
	picotrack->SetFlag(mcFlag);
	Short_t genIndex = mcpart->GetGeneratorIndex();
	picotrack->SetGeneratorIndex(genIndex);
      }
    }

    ++nacc;
  }
}

//________________________________________________________________________
AliVParticle* AliEmcalPicoTrackMaker::GetMCParticle(Int_t label) 
{
  if (!fMCParticles) return 0;
  Int_t index = label;
  if (fMCParticlesMap) index = fMCParticlesMap->At(label);
  if (index < 0 || index >= fMCParticles->GetEntriesFast()) return 0;
  AliVParticle *part = static_cast<AliVParticle*>(fMCParticles->At(index));
  return part;
}
