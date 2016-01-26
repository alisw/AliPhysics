//
// Task to filter Esd tracks and propagate to Emcal surface.
//
// Author: C.Loizides

#include "AliEmcalEsdTrackFilterTask.h"
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TGeoGlobalMagField.h>
#include <AliAnalysisManager.h>
#include <AliEMCALRecoUtils.h>
#include <AliESDEvent.h>
#include <AliESDtrackCuts.h>
#include <AliMagF.h>
#include <AliTrackerBase.h>


ClassImp(AliEmcalEsdTrackFilterTask)

//________________________________________________________________________
AliEmcalEsdTrackFilterTask::AliEmcalEsdTrackFilterTask() : 
  AliAnalysisTaskSE("AliEmcalEsdTrackFilterTask"),
  fEsdTrackCuts(0),
  fDoSpdVtxCon(0),
  fHybridTrackCuts(0),
  fTracksName(),
  fIncludeNoITS(kTRUE),
  fDoPropagation(kFALSE),
  fDist(440),
  fTrackEfficiency(0),
  fIsMC(kFALSE),
  fEsdEv(0),
  fTracks(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalEsdTrackFilterTask::AliEmcalEsdTrackFilterTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fEsdTrackCuts(0),
  fDoSpdVtxCon(0),
  fHybridTrackCuts(0),
  fTracksName("EsdTracksOut"),
  fIncludeNoITS(kTRUE),
  fDoPropagation(kFALSE),
  fDist(440),
  fTrackEfficiency(0),
  fIsMC(kFALSE),
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
AliEmcalEsdTrackFilterTask::~AliEmcalEsdTrackFilterTask()
{
  //Destructor

  delete fEsdTrackCuts;
}

//________________________________________________________________________
void AliEmcalEsdTrackFilterTask::UserCreateOutputObjects()
{
  // Create histograms.

  fTracks = new TClonesArray("AliESDtrack");
  fTracks->SetName(fTracksName);


  if (!fEsdTrackCuts) {
    if (fDoSpdVtxCon) {
      AliInfo("No track cuts given, creating default (standard only TPC) cuts");
      fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      fEsdTrackCuts->SetPtRange(0.15,1e3);
    } 
    else {
      AliWarning("No track cuts given, but maybe this is indeed intended?");
    }
  }
}

//________________________________________________________________________
void AliEmcalEsdTrackFilterTask::UserExec(Option_t *) 
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
  fTracks->Delete();
  if (!(InputEvent()->FindListObject(fTracksName)))
    InputEvent()->AddObject(fTracks);

  if (!fHybridTrackCuts) { // constrain TPC tracks to SPD vertex if fDoSpdVtxCon==kTRUE
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

        if (fTrackEfficiency) {
          Double_t r = gRandom->Rndm();
          if (fTrackEfficiency->Eval(ntrack->Pt()) < r)
            continue;
        }

        if (fDoPropagation)
          AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(ntrack,fDist);
        new ((*fTracks)[ntrnew++]) AliESDtrack(*ntrack);
        delete ntrack;
      }
    } else { /* no spd vtx constraint */
      Int_t ntr = fEsdEv->GetNumberOfTracks();
      for (Int_t i=0, ntrnew=0; i<ntr; ++i) {
        AliESDtrack *etrack = fEsdEv->GetTrack(i);
        if (!etrack)
          continue;

        if ((fEsdTrackCuts!=0) && !fEsdTrackCuts->AcceptTrack(etrack))
          continue;

        if (fTrackEfficiency) {
          Double_t r = gRandom->Rndm();
          if (fTrackEfficiency->Eval(etrack->Pt()) < r)
            continue;
        }

        AliESDtrack *ntrack = new ((*fTracks)[ntrnew++]) AliESDtrack(*etrack);
        if (fDoPropagation)
          AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(ntrack,fDist);
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
        if (fTrackEfficiency) {
          Double_t r = gRandom->Rndm();
          if (fTrackEfficiency->Eval(etrack->Pt()) < r)
            continue;
        }
        AliESDtrack *newTrack = new ((*fTracks)[ntrnew]) AliESDtrack(*etrack);
        if (fDoPropagation)
          AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(newTrack,fDist);
        newTrack->SetBit(BIT(22),0); 
        newTrack->SetBit(BIT(23),0);
        if (!fMCEvent) newTrack->SetLabel(0);
        ++ntrnew;
      } else if (fHybridTrackCuts->AcceptTrack(etrack)) {
        if (!etrack->GetConstrainedParam())
          continue;
        UInt_t status = etrack->GetStatus();
        if (!fIncludeNoITS && ((status&AliESDtrack::kITSrefit)==0))
          continue;

        if (fTrackEfficiency) {
          Double_t r = gRandom->Rndm();
          if (fTrackEfficiency->Eval(etrack->Pt()) < r)
            continue;
        }
        AliESDtrack *newTrack = new ((*fTracks)[ntrnew]) AliESDtrack(*etrack);
        if (fDoPropagation)
          AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(newTrack,fDist);
        const AliExternalTrackParam* constrainParam = etrack->GetConstrainedParam();
        newTrack->Set(constrainParam->GetX(),
            constrainParam->GetAlpha(),
            constrainParam->GetParameter(),
            constrainParam->GetCovariance());
        if ((status&AliESDtrack::kITSrefit)==0) {
          newTrack->SetBit(BIT(22),0); //type 2
          newTrack->SetBit(BIT(23),1);
        } else {
          newTrack->SetBit(BIT(22),1); //type 1
          newTrack->SetBit(BIT(23),0);
        }
        if (!fMCEvent) newTrack->SetLabel(0);
        ++ntrnew;
      }
    }
  }
}
