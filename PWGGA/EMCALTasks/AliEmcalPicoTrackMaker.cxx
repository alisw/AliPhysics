// $Id: AliEmcalPicoTrackMaker.cxx  $
//
// Class to make PicoTracks in AOD/ESD events.
//
//

#include <TClonesArray.h>

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliPicoTrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"

#include "AliEmcalPicoTrackMaker.h"

ClassImp(AliEmcalPicoTrackMaker)

//________________________________________________________________________
AliEmcalPicoTrackMaker::AliEmcalPicoTrackMaker() : 
  AliAnalysisTaskSE("AliEmcalPicoTrackMaker"),
  fESDtrackCuts(0),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalPicoTrackMaker::AliEmcalPicoTrackMaker(const char *name) : 
  AliAnalysisTaskSE("AliEmcalPicoTrackMaker"),
  fESDtrackCuts(0),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.

  if (!name)
    return;

  SetName(name);

  fBranchNames = "ESD:AliESDHeader.,AliESDRun.,SPDVertex.,Tracks";
}

//________________________________________________________________________
AliEmcalPicoTrackMaker::~AliEmcalPicoTrackMaker()
{
  //Destructor

}

//________________________________________________________________________
void AliEmcalPicoTrackMaker::UserCreateOutputObjects()
{
  // Create histograms.

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

  RetrieveEventObjects();

  // add tracks to event if not yet there
  if (!(InputEvent()->FindListObject(fTracksOutName)))
    InputEvent()->AddObject(fTracksOut);

  Int_t Ntracks = GetNumberOfTracks();
  for (Int_t iTracks = 0, nacc = 0; iTracks < Ntracks; ++iTracks) {
    AliVTrack *track = GetTrack(iTracks);
    
    if (!AcceptTrack(track))
      continue;

    Int_t label = 0;

    if (track->InheritsFrom("AliAODTrack")) {
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
      if (aodtrack->TestFilterBit(fAODfilterBits[0]))
	label = 0;
      else if (aodtrack->TestFilterBit(fAODfilterBits[1]))
	label = 1;
      else //if (aodtrack->TestFilterBit(fAODfilterBits[2]))
	label = 2;
    }
    else {
      label = track->GetLabel();
    }

    AliPicoTrack *picotrack = new ((*fTracksOut)[nacc]) AliPicoTrack(track->Pt(), track->Eta(), track->Phi(), 
					   track->Charge(), label, 
					   track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal(), track->IsEMCAL());
    if (track->IsEMCAL()) {
      picotrack->SetEMCALcluster(track->GetEMCALcluster());
    }
    //cout << iTracks << " - is emcal = " << track->IsEMCAL() << ", phiemc = " << track->GetTrackPhiOnEMCal() << ", etaemc = " << track->GetTrackEtaOnEMCal() << ", emcid = " << track->GetEMCALcluster() << endl;
    ++nacc;
  }
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackMaker::AcceptTrack(AliVTrack *track)
{
  if (!track)
    return kFALSE;

  if (track->InheritsFrom("AliAODTrack")) {
    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
    if (aodtrack) {
      //cout << "filter bit = " << fFilterBit << ", filter map = " << aodtrack->GetFilterMap() << endl;
      //return aodtrack->TestFilterBit(fAODfilterBits[0]+fAODfilterBits[1]+fAODfilterBits[2]);
      return aodtrack->IsHybridGlobalConstrainedGlobal();
    }
    else {
      AliError("Could not cast AOD track!");
      return kFALSE;
    }
  }
  else if (track->InheritsFrom("AliESDtrack")) { 
    if (fESDtrackCuts) {
      AliESDtrack *esdtrack = dynamic_cast<AliESDtrack*>(track);
      if (esdtrack) {
	return fESDtrackCuts->AcceptTrack(esdtrack);
      }
      else {
	AliError("Could not cast ESD track!");
	return kFALSE;
      }
    }
    else {
      return kTRUE;
    }
  }
  else if (track->InheritsFrom("PicoTrack")) {
    AliWarning("PicoTrack: nothing to filter!");
    return kTRUE;
  }
}

//________________________________________________________________________
void AliEmcalPicoTrackMaker::RetrieveEventObjects()
{
  fTracksIn = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksInName));

  if (!fTracksIn) {
    AliError(Form("ERROR: Could not retrieve tracks %s!", fTracksInName.Data())); 
  }
}

//________________________________________________________________________
AliVTrack* AliEmcalPicoTrackMaker::GetTrack(const Int_t i) const
{
  if (fTracksIn)
    return dynamic_cast<AliVTrack*>(fTracksIn->At(i));
  else
    return 0;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackMaker::GetNumberOfTracks() const
{
  if (fTracksIn)
    return fTracksIn->GetEntriesFast();
  else
    return 0;
}
