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
  fMinTrackPt(0),
  fMaxTrackPt(1000),
  fMinTrackEta(-0.9),
  fMaxTrackEta(0.9),
  fMinTrackPhi(-10),
  fMaxTrackPhi(10),
  fTrackEfficiency(1),
  fIncludeNoITS(kTRUE),
  fUseNegativeLabels(kFALSE),
  fIsMC(kFALSE),
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
  fESDtrackCuts(0),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fMinTrackPt(0),
  fMaxTrackPt(1000),
  fMinTrackEta(-0.9),
  fMaxTrackEta(0.9),
  fMinTrackPhi(-10),
  fMaxTrackPhi(10),
  fTrackEfficiency(1),
  fIncludeNoITS(kTRUE),
  fUseNegativeLabels(kFALSE),
  fIsMC(kFALSE),
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

    if (track->Pt() > fMaxTrackPt || track->Pt() < fMinTrackPt)
      continue;

    if (track->Eta() < fMinTrackEta || track->Eta() > fMaxTrackEta || 
	track->Phi() < fMinTrackPhi || track->Phi() > fMaxTrackPhi)
      continue;

    Bool_t isEmc = kFALSE;
    Int_t type = -1;
    if (esdMode) {
      AliESDtrack *esdtrack = static_cast<AliESDtrack*>(track);
      if (fESDtrackCuts && !fESDtrackCuts->AcceptTrack(esdtrack))
	continue;
      type = esdtrack->GetTRDNchamberdEdx();
      if (!fIncludeNoITS && (type==2))
	continue;
      isEmc = track->IsEMCAL();
    } else {
      AliAODTrack *aodtrack = static_cast<AliAODTrack*>(track);
      if (fAODfilterBits[0] < 0) {
	if (aodtrack->IsHybridGlobalConstrainedGlobal())
	  type = 3;
	else /*not a good track*/
	  continue;
      }
      else {
	if (aodtrack->TestFilterBit(fAODfilterBits[0])) {
	  type = 0;
	}
	else if (aodtrack->TestFilterBit(fAODfilterBits[1])) {
	  if ((aodtrack->GetStatus()&AliESDtrack::kITSrefit)==0) {
	    if (fIncludeNoITS)
	      type = 2;
	    else
	      continue;
	  }
	  else {
	    type = 1;
	  }
	}
	else {/*not a good track*/
	  continue;
	}
      }

      if (TMath::Abs(track->GetTrackEtaOnEMCal()) < 0.75 && 
	  track->GetTrackPhiOnEMCal() > 70 * TMath::DegToRad() && 
	  track->GetTrackPhiOnEMCal() < 190 * TMath::DegToRad())
	isEmc = kTRUE;
    }

    if (fTrackEfficiency < 1) {
      Double_t r = gRandom->Rndm();
      if (fTrackEfficiency < r) 
	continue;
    }

    Int_t label = 0;
    if (fIsMC) {
      if (track->GetLabel() > 0) {
	label = track->GetLabel();
      }
      else {
	if (!fUseNegativeLabels)
	  label = -track->GetLabel();
      }
      
      if (label == 0) 
	label = 99999;
    }

    /*AliPicoTrack *picotrack =*/ new ((*fTracksOut)[nacc]) AliPicoTrack(track->Pt(), 
									 track->Eta(), 
									 track->Phi(), 
									 track->Charge(), 
									 label,
									 type,
									 track->GetTrackEtaOnEMCal(), 
									 track->GetTrackPhiOnEMCal(), 
									 isEmc);
    ++nacc;
  }
}
