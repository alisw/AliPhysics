// $Id$
//
// Task to constrain TPC tracks to SPD vertex.
//
// Author: C.Loizides

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
  AliAnalysisTaskSE("TpcSpdVertexConstrainedTracks"),
  fEsdTrackCuts(0),
  fDoSpdVtxCon(0),
  fHybridTrackCuts(0),
  fTracksName(),
  fEsdEv(0),
  fTracks(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalEsdTpcTrackTask::AliEmcalEsdTpcTrackTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fEsdTrackCuts(0),
  fDoSpdVtxCon(0),
  fHybridTrackCuts(0),
  fTracksName("TpcSpdVertexConstrainedTracks"),
  fEsdEv(0),
  fTracks(0)
{
  // Constructor.

  if (!name)
    return;

  SetName(name);

  fBranchNames = "ESD:AliESDHeader.,AliESDRun.,SPDVertex.,Tracks";
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
  // Create histograms.

  fTracks = new TClonesArray("AliESDtrack");
  fTracks->SetName(fTracksName);

  if (!fEsdTrackCuts) {
    AliInfo("No track cuts given, creating default (standard only TPC) cuts");
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

  if (!fHybridTrackCuts) { // contrain TPC tracks to SPD vertex if fDoSpdVtxCon==kTRUE
    am->LoadBranch("AliESDRun.");
    am->LoadBranch("AliESDHeader.");
    am->LoadBranch("Tracks");

    if (fDoSpdVtxCon) {
      if (!TGeoGlobalMagField::Instance()->GetField()) { // construct field map
        fEsdEv->InitMagneticField();
      }
      am->LoadBranch("SPDVertex.");
      const AliESDVertex *vtxSPD = fEsdEv->GetPrimaryVertexSPD();
      if (!vtxSPD) {
        AliError("No SPD vertex, returning");
        return;
      }
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
    } else { /* no spd vtx constraint */
      Int_t ntr = fEsdEv->GetNumberOfTracks();
      for (Int_t i=0, ntrnew=0; i<ntr; ++i) {
        AliESDtrack *etrack = fEsdEv->GetTrack(i);
        if (!etrack)
          continue;
        if (!fEsdTrackCuts->AcceptTrack(etrack))
          continue;
        new ((*fTracks)[ntrnew++]) AliESDtrack(*etrack);
      }
    }

  } else { // use hybrid track cuts

    am->LoadBranch("Tracks");
    Int_t ntr = fEsdEv->GetNumberOfTracks();
    for (Int_t i=0, ntrnew=0; i<ntr; ++i) {
      AliESDtrack *etrack = fEsdEv->GetTrack(i);
      if (!etrack)
        continue;
      if (fEsdTrackCuts->AcceptTrack(etrack)) {
        new ((*fTracks)[ntrnew]) AliESDtrack(*etrack);
        AliESDtrack *newTrack = static_cast<AliESDtrack*>(fTracks->At(ntrnew));
        newTrack->SetLabel(0);
        ++ntrnew;
      } else if (fHybridTrackCuts->AcceptTrack(etrack)) {
        if (etrack->GetConstrainedParam()) {
          new ((*fTracks)[ntrnew]) AliESDtrack(*etrack);
          AliESDtrack *newTrack = static_cast<AliESDtrack*>(fTracks->At(ntrnew));
          ++ntrnew;
          const AliExternalTrackParam* constrainParam = etrack->GetConstrainedParam();
          newTrack->Set(constrainParam->GetX(),
                        constrainParam->GetAlpha(),
                        constrainParam->GetParameter(),
                        constrainParam->GetCovariance());
          UInt_t status = etrack->GetStatus();
          if ((status&AliESDtrack::kITSrefit)==0)
            newTrack->SetLabel(2);
          else
            newTrack->SetLabel(1);
        }
      }
    }
  }
}
