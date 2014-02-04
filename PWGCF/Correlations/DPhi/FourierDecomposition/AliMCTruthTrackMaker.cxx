//
// Select MC tracks, feed them via PicoTrack-like list

#include "AliMCTruthTrackMaker.h"

#include <TClonesArray.h>
#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliLog.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliESDEvent.h"
#include "AliStack.h"

ClassImp(AliMCTruthTrackMaker)

//________________________________________________________________________
AliMCTruthTrackMaker::AliMCTruthTrackMaker() :
AliAnalysisTaskSE("AliMCTruthTrackMaker"),
fTracksOutName("MCTruthTracks"),
fChargedMC(kFALSE),
fFillMuMothers(kFALSE),
fTriggerMatch(kTRUE),
fEtaMax(1.0),
fInit(kFALSE),
fEsdMode(kFALSE),
fTracksIn(0x0),
fTracksOut(0x0),
fESD(0x0),
fMC(0x0),
fStack(0x0)
{
  // Constructor.
}

//________________________________________________________________________
AliMCTruthTrackMaker::AliMCTruthTrackMaker(const char *name) :
AliAnalysisTaskSE(name),
fTracksOutName("MCTruthTracks"),
fChargedMC(kFALSE),
fFillMuMothers(kFALSE),
fTriggerMatch(kTRUE),
fEtaMax(1.0),
fInit(kFALSE),
fEsdMode(kFALSE),
fTracksIn(0x0),
fTracksOut(0x0),
fESD(0x0),
fMC(0x0),
fStack(0x0)
{
  // Constructor.
}

//________________________________________________________________________
AliMCTruthTrackMaker::~AliMCTruthTrackMaker()
{
  // Destructor.
}

//________________________________________________________________________
void AliMCTruthTrackMaker::UserCreateOutputObjects()
{
  // Create my user objects.
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliFatal("No analysis manager!");
    return;
  }
  
  AliVEventHandler *handler = mgr->GetInputEventHandler();
  if (!handler) {
    AliFatal("No input handler!");
    return;
  }
  
  if (handler->InheritsFrom("AliESDInputHandler"))
    fEsdMode = kTRUE;
  else
    fEsdMode = kFALSE;
  
  if (fEsdMode)
    fTracksOut = new TClonesArray("AliMCParticle");
  else
    fTracksOut = new TClonesArray("AliAODMCParticle");
  fTracksOut->SetName(fTracksOutName);
  
}

//________________________________________________________________________
void AliMCTruthTrackMaker::UserExec(Option_t *)
{
  // Main loop, called for each event.
  if (!InputEvent()) {
    AliError("Could not retrieve event! Returning");
    return;
  }
  
  if (fEsdMode && !MCEvent()) {
    AliError("Could not retrieve ESD MC event! Returning");
    return;
  }
  
  if (!fInit) {
    // add tracks to event if not yet there
    if (!(InputEvent()->FindListObject(fTracksOutName)))
      InputEvent()->AddObject(fTracksOut);
    
    if (!fEsdMode) {
      fTracksIn = static_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fTracksIn) {
        AliError("Could not retrieve AOD MC particles! Returning");
        return;
      }
      TClass *cl = fTracksIn->GetClass();
      if (!cl->GetBaseClass("AliAODMCParticle")) {
        AliError(Form("%s: Collection %s does not contain AliAODMCParticle!", GetName(), AliAODMCParticle::StdBranchName()));
        fTracksIn = 0;
        return;
      }
    }
    
    fInit = kTRUE;
  }
  
  // clear container (normally a null operation as the event should clean it already)
  fTracksOut->Delete();
  
  Int_t nacc = 0;
  
  // first, fill MC tracks up to fEtaMax
  const Int_t Ntracks = GetNumberOfTracks();
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
    AliVParticle* track = GetTrack(iTracks);

    if (!track)
      continue;
    
    if (TMath::Abs(track->Eta()) > fEtaMax)
      continue;

    if (fChargedMC && track->Charge() == 0) continue;

    AddTrack(track, nacc);
    
    ++nacc;
  }
  
  // now loop over reconstructed tracks and if accepted by muon cuts, add their primary mother
  if (fFillMuMothers) {
    if (!fEsdMode) {
      AliError("Muon mother particles can only be filled when running on ESD! Returning");
      return;
    }
    
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) { AliError("ESD event not found. Nothing done!"); return; }
    Int_t ntrks = fESD->GetNumberOfMuonTracks();

    fMC = MCEvent();
    if (!fMC) { AliError("MC event not available."); return; }
    fStack = fMC->Stack();
    
    for (Int_t iTrack = 0; iTrack<ntrks; iTrack++)
    {
      Int_t label = 0;
      AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iTrack);
      if(muonTrack)
      {
        if (IsGoodMUONtrack(*muonTrack)) {
          label =  TMath::Abs(muonTrack->GetLabel());
          if (label>=fMC->GetNumberOfTracks()) {
            AliError(Form("Label %d larger than number of particles on stack %d\n",label,fMC->GetNumberOfTracks()));
            continue;
          }
          AliVParticle * fpmMu = (AliVParticle *) fMC->GetTrack(GetFirstPrimaryMother(label));
          if (fpmMu) {
            AddTrack(fpmMu, nacc);
            ++nacc;
          } else {
            AliError("didn't receive a valid first primary mother here");
          }
        }
      }
    }
    
  }
}

//________________________________________________________________________
Int_t AliMCTruthTrackMaker::GetNumberOfTracks() const
{
  if (fEsdMode)
    return MCEvent()->GetNumberOfTracks();
  else
    return fTracksIn->GetEntries();
}

//________________________________________________________________________
AliVParticle* AliMCTruthTrackMaker::GetTrack(Int_t i)
{
  if (fEsdMode) {
    if (!MCEvent()->IsPhysicalPrimary(i))
      return 0;
    
    return MCEvent()->GetTrack(i);
  }
  else {
    AliAODMCParticle *part = static_cast<AliAODMCParticle*>(fTracksIn->At(i));
    if (!part->IsPhysicalPrimary())
      return 0;
    
    return part;
  }
}

//________________________________________________________________________
void AliMCTruthTrackMaker::AddTrack(AliVParticle *track, Int_t nacc)
{
  if (fEsdMode) {
    AliMCParticle *part = static_cast<AliMCParticle*>(track);
    new ((*fTracksOut)[nacc]) AliMCParticle(part->Particle(), 0, part->Label());
  }
  else {
    AliAODMCParticle *part = static_cast<AliAODMCParticle*>(track);
    new ((*fTracksOut)[nacc]) AliAODMCParticle(*part);
  }
}

//________________________________________________________________________
Bool_t AliMCTruthTrackMaker::IsGoodMUONtrack(AliESDMuonTrack &track)
{
  // Applying track cuts for MUON tracks
  
  if (!track.ContainTrackerData())
    return kFALSE;
  Double_t thetaTrackAbsEnd = TMath::ATan(track.GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
  if ((thetaTrackAbsEnd < 2.) || (thetaTrackAbsEnd > 10.))
    return kFALSE;
  Double_t eta = track.Eta();
  if ((eta < -4.) || (eta > -2.5))
    return kFALSE;
  if (fTriggerMatch) {
    if (!track.ContainTriggerData())
      return kFALSE;
    if (track.GetMatchTrigger() < 0.5)
      return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMCTruthTrackMaker::IsGoodMUONtrack(AliAODTrack &track)
{
  // Applying track cuts for MUON tracks
  
  if (!track.IsMuonTrack())
    return kFALSE;
  Double_t dThetaAbs = TMath::ATan(track.GetRAtAbsorberEnd()/505.)* TMath::RadToDeg();
  if ((dThetaAbs<2.) || (dThetaAbs>10.))
    return kFALSE;
  Double_t dEta = track.Eta();
  if ((dEta<-4.) || (dEta>-2.5))
    return kFALSE;
  if (fTriggerMatch) {
    if (track.GetMatchTrigger()<0.5)
      return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Int_t AliMCTruthTrackMaker::GetFirstPrimaryMother(Int_t muonlabel)
{
  AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(muonlabel);
  if(McParticle->GetMother()<fStack->GetNprimary()) return McParticle->GetMother();
  else
  {
    Int_t motherlabel = McParticle->GetMother();
    while(motherlabel > -1)
    {
      AliMCParticle *MotherParticle  = (AliMCParticle*)fMC->GetTrack(motherlabel);
      if(MotherParticle->GetMother()<fStack->GetNprimary()) break;
      else motherlabel = MotherParticle->GetMother();
    }
    AliMCParticle *FirstSecondaryMotherParticle  = (AliMCParticle*)fMC->GetTrack(motherlabel);
    return FirstSecondaryMotherParticle->GetMother();
  }
}
