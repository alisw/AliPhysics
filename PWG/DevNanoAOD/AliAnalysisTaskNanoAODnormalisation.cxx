#include "AliAnalysisTaskNanoAODnormalisation.h"
#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>
#include <AliESDtrackCuts.h>
#include <AliInputEventHandler.h>
#include <AliNanoFilterNormalisation.h>
#include <AliVEvent.h>

#include <TChain.h>
#include <TObject.h>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {
  constexpr char kNames[2][16 ]{"Filter","skimming"};
}

AliAnalysisTaskNanoAODnormalisation::AliAnalysisTaskNanoAODnormalisation(std::string taskName) : AliAnalysisTaskSE{taskName.data()},
  fFirst{true},
  fOutputList{nullptr},
  fCandidateEvents{nullptr},
  fSelectedEvents{nullptr}
{
  DefineInput(0, TChain::Class());  
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskNanoAODnormalisation::~AliAnalysisTaskNanoAODnormalisation() {
  if (fOutputList)
    delete fOutputList;
}

void AliAnalysisTaskNanoAODnormalisation::UserCreateOutputObjects() {
  fOutputList = new TList;
  fOutputList->SetOwner(true);
  int nMultBins = 100;
  float multBegin = 0;
  float multEnd = 100;

  for (int iF{0}; iF < 2; ++iF) {
    fCandidateEvents[iF] = new TH2D(Form("fCandidateEvents%s",kNames[iF]), ";Multiplicity estimator;", nMultBins, multBegin, multEnd, 5, -0.5, 4.5);
    fSelectedEvents[iF]  = (TH2D*) fCandidateEvents[iF]->Clone(Form("fSelectedEvents%s",kNames[iF]));

    std::string labels[5]{"Input events", "Triggered events", "Triggered + Quality cuts", "Triggered + QC + Reco vertex", "Analysis events"};
    for (int iType{0}; iType < 5; ++iType) {
      fCandidateEvents[iF]->GetYaxis()->SetBinLabel(iType + 1, labels[iType].data());
      fSelectedEvents[iF]->GetYaxis()->SetBinLabel(iType + 1, labels[iType].data());
    }
    fOutputList->Add(fCandidateEvents[iF]);
    fOutputList->Add(fSelectedEvents[iF]);
  }
  PostData(1, fOutputList);
}

Bool_t AliAnalysisTaskNanoAODnormalisation::UserNotify() {
  AliAODInputHandler* handler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!handler) {
    ::Warning("AliAnalysisTaskNanoAODnormalisation::UserNotify","Missing input handler");
    return true;
  }
  TList* userInfo = handler->GetUserInfo();
  if (!userInfo) {
    ::Warning("AliAnalysisTaskNanoAODnormalisation::UserNotify","Missing User Info");
    return true;
  }

  for (int iF{0}; iF < 2; ++iF) {
    AliNanoFilterNormalisation* normalisation = (AliNanoFilterNormalisation*)userInfo->FindObject(Form("NanoAOD%s_scaler",kNames[iF]));
    if (!normalisation) {
      if (!iF)
        ::Fatal("AliAnalysisTaskNanoAODnormalisation::UserNotify","Missing filtering scalers!");
      else
        ::Warning("AliAnalysisTaskNanoAODnormalisation::UserNotify","Missing skimming scalers!");
      return true;
    }
    if (fFirst) {
      /// We do the copy, in case someone wanted a different multiplicity binning
      normalisation->GetCandidateEventsHistogram()->Copy(*fCandidateEvents[iF]);
      normalisation->GetSelectedEventsHistogram()->Copy(*fSelectedEvents[iF]);
      fCandidateEvents[iF]->SetName(Form("fCandidateEvents_%s",kNames[iF]));
      fSelectedEvents[iF]->SetName(Form("fSelectedEvents_%s",kNames[iF]));
      if (iF)
        fFirst = false;
    } else {
      fCandidateEvents[iF]->Add(normalisation->GetCandidateEventsHistogram());
      fSelectedEvents[iF]->Add(normalisation->GetSelectedEventsHistogram());
    }
  }
  PostData(1, fOutputList);  
  return true;
}

AliAnalysisTaskNanoAODnormalisation* AliAnalysisTaskNanoAODnormalisation::AddTask(std::string name) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisTaskNanoAODnormalisation::AddTask", "No analysis manager to connect to.");
    return nullptr;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AliAnalysisTaskNanoAODnormalisation::AddTask", "This task requires an input event handler");
    return nullptr;
  }

  AliAnalysisTaskNanoAODnormalisation *task = new AliAnalysisTaskNanoAODnormalisation(name);
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Skimming_Normalisation", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(), "NanoAODskimmingNormalisation"));
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
