// $Id$
//
// Class to select tracks in MC events.
//
// Author: S. Aiola

#include "AliEmcalMCTrackSelector.h"

#include <TClonesArray.h>

#include "AliNamedArrayI.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliLog.h"

ClassImp(AliEmcalMCTrackSelector)

//________________________________________________________________________
AliEmcalMCTrackSelector::AliEmcalMCTrackSelector() : 
  AliAnalysisTaskSE("AliEmcalMCTrackSelector"),
  fTracksOutName("PicoTracks"),
  fRejectNK(kFALSE),
  fChargedMC(kFALSE),
  fTracksMapName(""),
  fEtaMax(1),
  fInit(kFALSE),
  fEsdMode(kFALSE),
  fTracksIn(0),
  fTracksOut(0),
  fTracksMap(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalMCTrackSelector::AliEmcalMCTrackSelector(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksOutName("PicoTracks"),
  fRejectNK(kFALSE),
  fChargedMC(kFALSE),
  fTracksMapName(""),
  fEtaMax(1),
  fInit(kFALSE),
  fEsdMode(kFALSE),
  fTracksIn(0),
  fTracksOut(0),
  fTracksMap(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalMCTrackSelector::~AliEmcalMCTrackSelector()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalMCTrackSelector::UserCreateOutputObjects()
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

  fTracksMapName = fTracksOutName;
  fTracksMapName += "_Map";
  fTracksMap = new AliNamedArrayI(fTracksMapName, 99999);
}

//________________________________________________________________________
void AliEmcalMCTrackSelector::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  if (!InputEvent()) {
    AliError("Could not retrieve event! Returning");
    return;
  }

  if (fEsdMode && !MCEvent()) {
    AliError("Could not retrieve MC event! Returning");
    return;
  }

  if (!fInit) {
    // add tracks to event if not yet there
    if (!(InputEvent()->FindListObject(fTracksOutName))) 
      InputEvent()->AddObject(fTracksOut);

    if (!(InputEvent()->FindListObject(fTracksMapName)))
      InputEvent()->AddObject(fTracksMap);

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
  fTracksMap->Clear();

  const Int_t Ntracks = GetNumberOfTracks();

  // loop over tracks
  for (Int_t iTracks = 0, nacc = 0; iTracks < Ntracks; ++iTracks) {

    if (iTracks >= fTracksMap->GetSize())
      fTracksMap->Set(iTracks*2);

    fTracksMap->AddAt(-1, iTracks);

    AliVParticle* track = GetTrack(iTracks);

    if (!track)
      continue;

    if (TMath::Abs(track->Eta()) > fEtaMax) 
      continue;
    
    Int_t pdgCode = track->PdgCode();
    if (fRejectNK && (pdgCode == 130 || pdgCode == 2112)) continue;
    
    if (fChargedMC && track->Charge() == 0) continue;

    fTracksMap->AddAt(nacc, iTracks);

    AddTrack(track, nacc);

    ++nacc;
  }
}

//________________________________________________________________________
Int_t AliEmcalMCTrackSelector::GetNumberOfTracks() const
{
  if (fEsdMode)
    return MCEvent()->GetNumberOfTracks();
  else
    return fTracksIn->GetEntries();
}

//________________________________________________________________________
AliVParticle* AliEmcalMCTrackSelector::GetTrack(Int_t i)
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
void AliEmcalMCTrackSelector::AddTrack(AliVParticle *track, Int_t nacc)
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
