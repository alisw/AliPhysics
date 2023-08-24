/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <THashList.h>
#include <THistManager.h>
#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCountITStracks.h"
#include "AliAnalysisUtils.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskCountITStracks)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskCountITStracks::AliAnalysisTaskCountITStracks() :
    AliAnalysisTaskSE(),
    fAnalysisUtils(nullptr),
    fHistos(nullptr)
{

}

AliAnalysisTaskCountITStracks::AliAnalysisTaskCountITStracks(const char *name) :
    AliAnalysisTaskSE(name),
    fAnalysisUtils(nullptr),
    fHistos(nullptr)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskCountITStracks::~AliAnalysisTaskCountITStracks() {
  if(fAnalysisUtils) delete fAnalysisUtils;
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskCountITStracks::UserCreateOutputObjects(){
  fAnalysisUtils = new AliAnalysisUtils;

  fHistos = new THistManager("CountITStracks");
  fHistos->CreateTH1("countITStracksACside", "Total ITS track count for A- and C-side", 2, -0.5, 1.5);
  fHistos->CreateTH1("countITStracksEventAside", "Multiplicity distribution ITS tracks in EMCAL area A-side", 101, 0.5, 100.5);
  fHistos->CreateTH1("countITStracksEventCside", "Multiplicity distribution ITS tracks in EMCAL area C-side", 101, 0.5, 100.5);
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskCountITStracks::UserExec(Option_t *){
  // idea: count number of ITS tracks for A- and C- side separately
  // (above pt-threshold) in the EMCAL area
  // In case of matching issues we should see an increase in the amount
  // of ITS stand-alone tracks.

  const int kAside = 0, kCside = 1;

  // Select events
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return;
  double vz = fInputEvent->GetPrimaryVertex()->GetZ();
  if(TMath::Abs(vz) > 10.) return;
  if(!fAnalysisUtils->IsVertexSelected2013pA(fInputEvent)) return;
  if(fAnalysisUtils->IsPileUpEvent(fInputEvent)) return;

  int nTracksAside(0), nTracksCside(0);
  int side(-1);
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); itrk++){
    AliVTrack *trk = static_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));

    // require track to be within the EMCAL area
    double eta = trk->Eta(), phi = trk->Phi();
    if(TMath::Abs(eta) > 0.5) continue;
    if(phi < 1.4 || phi > 3.1) continue;

    // check whether track is an ITS stand-alone track
    AliESDtrack *esdtrack(nullptr);
    AliAODTrack *aodtrack(nullptr);
    float dr, dz;
    if((esdtrack = static_cast<AliESDtrack *>(trk))) {
      if(!(esdtrack->GetStatus() & AliESDtrack::kITSpureSA)) continue;
      esdtrack->GetImpactParameters(dr, dz);
    } else if((aodtrack = static_cast<AliAODTrack *>(trk))) {
      if(!aodtrack->IsPureITSStandalone()) continue;
      aodtrack->GetImpactParameters(dr, dz);
    } else {
      continue;
    }

    // Check pointing to the primary vertex
    if(TMath::Abs(dr) > 3.2) continue;
    if(TMath::Abs(dz) > 2.) continue;

    // Apply pt-cut in order to reject low-pt ITS stand-alone tracks
    if(TMath::Abs(trk->Pt()) < 0.5) continue;

    if(eta < 0){
      side = kCside;
      nTracksCside++;
    }
    else{
      side = kAside;
      nTracksAside++;
    }

    fHistos->FillTH1("countITStracksACside", side);
  }

  // Multiplicities
  fHistos->FillTH1("countITStracksEventAside", nTracksAside);
  fHistos->FillTH1("countITStracksEventCside", nTracksCside);

  PostData(1, fHistos->GetListOfHistograms());
}

AliAnalysisTaskCountITStracks *AliAnalysisTaskCountITStracks::AddTaskCountITStracks(const char *name) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisTaskCountITStracks *task = new AliAnalysisTaskCountITStracks(name);
  mgr->AddTask(task);

  TString outputcont = mgr->GetCommonFileName();
  outputcont += ":CountITStracks";

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("histosCountITStracks", TList::Class(), AliAnalysisManager::kOutputContainer, outputcont));
  return task;
}
