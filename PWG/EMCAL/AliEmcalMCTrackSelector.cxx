// $Id$
//
// Class to select tracks in MC events.
//
// Author: S. Aiola

#include "AliEmcalMCTrackSelector.h"

#include <TClonesArray.h>
#include <TH1I.h>

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
  fInit(kFALSE),
  fTracksMapName(""),
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
  fInit(kFALSE),
  fTracksMapName(""),
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

  fTracksOut = new TClonesArray("AliMCParticle");
  fTracksOut->SetName(fTracksOutName);

  fTracksMapName = fTracksOutName;
  fTracksMapName += "_Map";
  fTracksMap = new TH1I(fTracksMapName, fTracksMapName, 1000, 0, 1);
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

  if (!fInit) {
    // add tracks to event if not yet there
    if (!(event->FindListObject(fTracksOutName))) 
      event->AddObject(fTracksOut);

    if (!(event->FindListObject(fTracksMapName)))
      event->AddObject(fTracksMap);

    fInit = kTRUE;
  }

  // clear container (normally a null operation as the event should clean it already)
  fTracksOut->Delete();

  const Int_t Ntracks = mcevent->GetNumberOfTracks();
  new (fTracksMap) TH1I(fTracksMapName, fTracksMapName, Ntracks-2, 0, 1);  // Ntracks - 2, we use also over- and uner-flow bins

  // loop over tracks
  for (Int_t iTracks = 0, nacc = 0; iTracks < Ntracks; ++iTracks) {

    fTracksMap->SetBinContent(iTracks, -1);

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

    fTracksMap->SetBinContent(iTracks, nacc);
    new ((*fTracksOut)[nacc]) AliMCParticle(track->Particle(), 0, track->Label());

    ++nacc;
  }
}
