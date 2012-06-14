// $Id$
//
// Class to make PicoTracks in AOD/ESD events.
//
// Author: S.Aiola, C.Loizides

#include <TClonesArray.h>
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVTrack.h"
#include "AliEmcalPicoTrackMaker.h"

ClassImp(AliEmcalPicoTrackMaker)

//________________________________________________________________________
AliEmcalPicoTrackMaker::AliEmcalPicoTrackMaker() : 
  AliAnalysisTaskSE("AliEmcalPicoTrackMaker"),
  fESDtrackCuts(0),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fMaxTrackPt(100),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.

  fAODfilterBits[0] = 0;
  fAODfilterBits[1] = 0;
}

//________________________________________________________________________
AliEmcalPicoTrackMaker::AliEmcalPicoTrackMaker(const char *name) : 
  AliAnalysisTaskSE(name),
  fESDtrackCuts(0),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fMaxTrackPt(100),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.

  fAODfilterBits[0] = 0;
  fAODfilterBits[1] = 0;
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
  if (!(InputEvent()->FindListObject(fTracksOutName))) {
    InputEvent()->AddObject(fTracksOut);
  }

  // clear container (normally a null operation as the event should clean it already)
  fTracksOut->Clear();

  // test if we are in ESD or AOD mode
  Bool_t esdMode = kTRUE;
  if (dynamic_cast<AliAODEvent*>(InputEvent())!=0)
    esdMode = kFALSE;
 
  // loop over tracks
  const Int_t Ntracks = fTracksIn->GetEntriesFast();
  for (Int_t iTracks = 0, nacc = 0; iTracks < Ntracks; ++iTracks) {

    AliVTrack *track = static_cast<AliVTrack*>(fTracksIn->At(iTracks));

    if (!track)
      continue;

    if (track->Pt() > fMaxTrackPt)
      continue;
    Bool_t isEmc = kFALSE;
    Int_t label = -1;
    if (esdMode) {
      if (fESDtrackCuts) {
        AliESDtrack *esdtrack = static_cast<AliESDtrack*>(track);
        if (!fESDtrackCuts->AcceptTrack(esdtrack))
          continue;
      }
      label = track->GetLabel();
      isEmc = track->IsEMCAL();
    } else {
      AliAODTrack *aodtrack = static_cast<AliAODTrack*>(track);
      if (aodtrack->TestFilterBit(fAODfilterBits[0]))
	label = 0;
      else if (aodtrack->TestFilterBit(fAODfilterBits[1]))
	label = 3;
      else /*not a good track*/
        continue;

      if (TMath::Abs(track->GetTrackEtaOnEMCal()) < 0.75 && 
	  track->GetTrackPhiOnEMCal() > 70 * TMath::DegToRad() && 
	  track->GetTrackPhiOnEMCal() < 190 * TMath::DegToRad())
	isEmc = kTRUE;
    }

    
    /*AliPicoTrack *picotrack =*/ new ((*fTracksOut)[nacc]) AliPicoTrack(track->Pt(), 
                                                                     track->Eta(), 
                                                                     track->Phi(), 
                                                                     track->Charge(), 
                                                                     label, 
                                                                     track->GetTrackEtaOnEMCal(), 
                                                                     track->GetTrackPhiOnEMCal(), 
                                                                     isEmc);
    ++nacc;
  }
}
