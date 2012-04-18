// $Id: AliEmcalAodTrackFilterTask.cxx  $
//
// Class to filter hybrid tracks in AOD events.
//
//

#include <TClonesArray.h>

#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliPicoTrack.h"

#include "AliEmcalAodTrackFilterTask.h"

ClassImp(AliEmcalAodTrackFilterTask)

//________________________________________________________________________
AliEmcalAodTrackFilterTask::AliEmcalAodTrackFilterTask() : 
  AliAnalysisTaskSE("AliEmcalAodTrackFilterTask"),
  fAODfilterBit(16),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fAOD(0),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalAodTrackFilterTask::AliEmcalAodTrackFilterTask(const char *name) : 
  AliAnalysisTaskSE("AliEmcalAodTrackFilterTask"),
  fAODfilterBit(16),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fAOD(0),
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
AliEmcalAodTrackFilterTask::~AliEmcalAodTrackFilterTask()
{
  //Destructor

}

//________________________________________________________________________
void AliEmcalAodTrackFilterTask::UserCreateOutputObjects()
{
  // Create histograms.

  fTracksOut = new TClonesArray("AliPicoTrack");
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

  RetrieveEventObjects();

  // add tracks to event if not yet there
  if (!(InputEvent()->FindListObject(fTracksOutName)))
    InputEvent()->AddObject(fTracksOut);

  Int_t Ntracks = GetNumberOfTracks();
  for (Int_t iTracks = 0, nacc = 0; iTracks < Ntracks; ++iTracks) {
    AliAODTrack *track = GetTrack(iTracks);

    if (!track)
      continue;
    
    if (!AcceptTrack(track))
      continue;

    AliPicoTrack *picotrack = new ((*fTracksOut)[nacc]) AliPicoTrack(track->Pt(), track->Eta(), track->Phi(), 
					   track->Charge(), track->GetLabel(), 
					   track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal(), track->IsEMCAL());
    if (track->IsEMCAL()) {
      picotrack->SetEMCALcluster(track->GetEMCALcluster());
    }
    //cout << iTracks << " - is emcal = " << track->IsEMCAL() << ", phiemc = " << track->GetTrackPhiOnEMCal() << ", etaemc = " << track->GetTrackEtaOnEMCal() << ", emcid = " << track->GetEMCALcluster() << endl;
    ++nacc;
  }
}


//________________________________________________________________________
Bool_t AliEmcalAodTrackFilterTask::AcceptTrack(AliAODTrack *track)
{
  AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
  if (aodtrack) {
    //cout << "filter bit = " << fFilterBit << ", filter map = " << aodtrack->GetFilterMap() << endl;
    return aodtrack->TestFilterBit(fAODfilterBit);
    //return aodtrack->IsHybridGlobalConstrainedGlobal();
  }
  return 1;
}

//________________________________________________________________________
void AliEmcalAodTrackFilterTask::SetRunPeriod(const char *p)
{
  if (!strcmp(p, "LHC11h")) {
    SetAODfilterBit(256+512); // hybrid tracks for LHC11h
  }
  else {
    AliWarning(Form("Run period %s not known. AOD filter bit not set.", p));
  }
}

//________________________________________________________________________
void AliEmcalAodTrackFilterTask::RetrieveEventObjects()
{
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    AliError("Task works only on AOD events, returning");
    return;
  }

  fTracksIn = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksInName));

  if (!fTracksIn) {
    AliError(Form("ERROR: Could not retrieve tracks %s!", fTracksInName.Data())); 
  }
}

//________________________________________________________________________
AliAODTrack* AliEmcalAodTrackFilterTask::GetTrack(const Int_t i) const
{
  return dynamic_cast<AliAODTrack*>(fTracksIn->At(i));
}

//________________________________________________________________________
Int_t AliEmcalAodTrackFilterTask::GetNumberOfTracks() const
{
  return fTracksIn->GetEntriesFast();
}
