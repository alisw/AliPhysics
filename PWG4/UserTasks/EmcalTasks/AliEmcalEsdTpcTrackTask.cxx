// $Id$

#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalEsdTpcTrackTask.h"
#include "AliMagF.h"

ClassImp(AliEmcalEsdTpcTrackTask)

//________________________________________________________________________
AliEmcalEsdTpcTrackTask::AliEmcalEsdTpcTrackTask() : 
  AliAnalysisTaskSE(),
  fEsdTrackCuts(0),
  fTracksName(),
  fEsdEv(0),
  fTracks(0)
{
  // Constructor
}

//________________________________________________________________________
AliEmcalEsdTpcTrackTask::AliEmcalEsdTpcTrackTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fEsdTrackCuts(0),
  fTracksName("TpcSpdVertexConstrainedTracks"),
  fEsdEv(0),
  fTracks(0)
{
  // Constructor
}

//________________________________________________________________________
AliEmcalEsdTpcTrackTask::~AliEmcalEsdTpcTrackTask()
{
  //Destructor

  delete fEsdTrackCuts;
}

//________________________________________________________________________
void AliEmcalEsdTpcTrackTask::UserCreateOutputObjects()
{
  // Create histograms

  fTracks = new TClonesArray("AliESDtrack");
  fTracks->SetName(fTracksName);

  if (!fEsdTrackCuts) {
    AliInfo("No track cuts given, creating default cuts");
    fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdTrackCuts->SetPtRange(0.15,1e3);
  }
}

//________________________________________________________________________
void AliEmcalEsdTpcTrackTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  fEsdEv = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fEsdEv) {
    AliError("Task works only on ESD events, returning");
    return;
  }

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }

  // add tracks to event if not yet there
  if (!(InputEvent()->FindListObject(fTracksName)))
    InputEvent()->AddObject(fTracks);

  am->LoadBranch("AliESDRun.");
  am->LoadBranch("AliESDHeader.");
  if (!TGeoGlobalMagField::Instance()->GetField()) { // construct field map
    const AliESDRun *erun = fEsdEv->GetESDRun();
    AliMagF *field = AliMagF::CreateFieldMap(erun->GetCurrentL3(),
                                             erun->GetCurrentDip(),
                                             AliMagF::kConvLHC,
                                             kFALSE,
                                             erun->GetBeamEnergy(),
                                             erun->GetBeamType());
    TGeoGlobalMagField::Instance()->SetField(field);
  }

  am->LoadBranch("SPDVertex.");
  const AliESDVertex *vtxSPD = fEsdEv->GetPrimaryVertexSPD();
  if (!vtxSPD) {
    AliError("No SPD vertex, returning");
    return;
  }

  am->LoadBranch("Tracks");
  Int_t ntr = fEsdEv->GetNumberOfTracks();
  for (Int_t i=0, ntrnew=0; i<ntr; ++i) {
    AliESDtrack *etrack = fEsdEv->GetTrack(i);
    if (!etrack)
      continue;
    if (!fEsdTrackCuts->AcceptTrack(etrack))
      continue;

    AliESDtrack *ntrack = AliESDtrackCuts::GetTPCOnlyTrack(fEsdEv,etrack->GetID());
    if (!ntrack)
      continue;
    if (ntrack->Pt()<=0) {
      delete ntrack;
      continue;
    }

    Double_t bfield[3] = {0,0,0};
    ntrack->GetBxByBz(bfield);
    AliExternalTrackParam exParam;
    Bool_t relate = ntrack->RelateToVertexBxByBz(vtxSPD,bfield,kVeryBig,&exParam);
    if (!relate) {
      delete ntrack;
      continue;
    }

    // set the constraint parameters to the track
    ntrack->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());

    if (ntrack->Pt()<=0) {
      delete ntrack;
      continue;
    }

    new ((*fTracks)[ntrnew++]) AliESDtrack(*ntrack);
    delete ntrack;
  }
}
