// $Id: AliEmcalAodTrackFilterTask.cxx | Fri Dec 6 10:29:20 2013 +0100 | Constantin Loizides  $
//
// Class to filter Aod tracks
//
// Author: C.Loizides

#include "AliEmcalAodTrackFilterTask.h"
#include <TClonesArray.h>
#include <TRandom3.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliAnalysisManager.h>
#include <AliEMCALRecoUtils.h>
#include <AliLog.h>

ClassImp(AliEmcalAodTrackFilterTask)

//________________________________________________________________________
AliEmcalAodTrackFilterTask::AliEmcalAodTrackFilterTask() : 
  AliAnalysisTaskSE("AliEmcalAodTrackFilterTask"),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fIncludeNoITS(kTRUE),
  fCutMaxFrShTPCClus(0.4),
  fUseNegativeLabels(kTRUE),
  fIsMC(kFALSE),
  fDoPropagation(kFALSE),
  fDist(440),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.

  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
}

//________________________________________________________________________
AliEmcalAodTrackFilterTask::AliEmcalAodTrackFilterTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fIncludeNoITS(kTRUE),
  fCutMaxFrShTPCClus(0.4),
  fUseNegativeLabels(kTRUE),
  fIsMC(kFALSE),
  fDoPropagation(kFALSE),
  fDist(440),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.

  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
  fBranchNames = "AOD:tracks";
}

//________________________________________________________________________
AliEmcalAodTrackFilterTask::~AliEmcalAodTrackFilterTask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalAodTrackFilterTask::UserCreateOutputObjects()
{
  // Create my user objects.

  fTracksOut = new TClonesArray("AliAODTrack");
  fTracksOut->SetName(fTracksOutName);
}

//________________________________________________________________________
void AliEmcalAodTrackFilterTask::UserExec(Option_t *) 
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

    AliAODTrack *track = static_cast<AliAODTrack*>(fTracksIn->At(iTracks));

    if (!track)
      continue;

    Int_t type = -1;
    if (fAODfilterBits[0] < 0) {
      if (track->IsHybridGlobalConstrainedGlobal())
	type = 3;
      else /*not a good track*/
	continue;
    } else {
      if (track->TestFilterBit(fAODfilterBits[0])) {
	type = 0;
      } else if (track->TestFilterBit(fAODfilterBits[1])) {
	if ((track->GetStatus()&AliVTrack::kITSrefit)==0) {
	  if (fIncludeNoITS)
	    type = 2;
	  else
	    continue;
	} else {
	  type = 1;
	}
      }
      else {/*not a good track*/
	continue;
      }
    }

    if (fCutMaxFrShTPCClus > 0) {
      Double_t frac = Double_t(track->GetTPCnclsS()) / Double_t(track->GetTPCncls());
      if (frac > fCutMaxFrShTPCClus) 
	continue;
    }

    AliAODTrack *newt = new ((*fTracksOut)[nacc]) AliAODTrack(*track);
    if (fDoPropagation)
      AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(newt,fDist);
    Int_t label = 0;
    if (fIsMC) {
      if (fUseNegativeLabels)
	label = track->GetLabel();
      else 
	label = TMath::Abs(track->GetLabel());
      if (label == 0) 
	AliDebug(2,Form("Track %d with label==0", iTracks));
    }
    if (type==0) {
      newt->SetBit(BIT(20),0);
      newt->SetBit(BIT(21),0);
    } else if (type==1) {
      newt->SetBit(BIT(20),1);
      newt->SetBit(BIT(21),0);
    } else if (type==2) {
      newt->SetBit(BIT(20),0);
      newt->SetBit(BIT(21),1);
    } else if (type==3) {
      newt->SetBit(BIT(20),1);
      newt->SetBit(BIT(21),1);
    }
    ++nacc;
  }
}
