// $Id$
//
// Class to select tracks in MC events.
//
// Author: S. Aiola

#include "AliEmcalMCTrackSelector.h"

#include <TClonesArray.h>

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliLog.h"

ClassImp(AliEmcalMCTrackSelector)

//________________________________________________________________________
AliEmcalMCTrackSelector::AliEmcalMCTrackSelector() : 
  AliAnalysisTaskSE("AliEmcalMCTrackSelector"),
  fTracksOutName("PicoTracks"),
  fRejectNK(kFALSE),
  fChargedMC(kFALSE),
  fTracksOut(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalMCTrackSelector::AliEmcalMCTrackSelector(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksOutName("PicoTracks"),
  fRejectNK(kFALSE),
  fChargedMC(kFALSE),
  fTracksOut(0)
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

  fTracksOut = new TClonesArray("AliMCParticle");
  fTracksOut->SetName(fTracksOutName);
}

//________________________________________________________________________
void AliEmcalMCTrackSelector::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  AliMCEvent *mcevent = MCEvent();
  AliVEvent *event = InputEvent();

  if (!event) {
    AliError("Could not retrieve event! Returning");
    return;
  }

  if (!mcevent) {
    AliError("Could not retrieve MC event! Returning");
    return;
  }

  // add tracks to event if not yet there
  if (!(event->FindListObject(fTracksOutName))) {
    event->AddObject(fTracksOut);
  }

  // clear container (normally a null operation as the event should clean it already)
  fTracksOut->Clear();
 
  // loop over tracks
  const Int_t Ntracks = mcevent->GetNumberOfTracks();
  for (Int_t iTracks = 0, nacc = 0; iTracks < Ntracks; ++iTracks) {

    if (!mcevent->IsPhysicalPrimary(iTracks))
      continue;

    AliMCParticle* track = dynamic_cast<AliMCParticle*>(mcevent->GetTrack(iTracks));

    if (!track)
      continue;

    if (TMath::Abs(track->Eta()) > 1) 
      continue;
    
    Int_t pdgCode = track->PdgCode();
    if (fRejectNK && (pdgCode == 130 || pdgCode == 2112)) continue;
    
    if (fChargedMC && track->Charge() == 0) continue;
    
    new ((*fTracksOut)[nacc]) AliMCParticle(track->Particle(), 0, track->Label());
    ++nacc;
  }
}
