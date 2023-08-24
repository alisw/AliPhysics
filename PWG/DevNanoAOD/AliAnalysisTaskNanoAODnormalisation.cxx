#include "AliAnalysisTaskNanoAODnormalisation.h"
#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>
#include <AliESDtrackCuts.h>
#include <AliInputEventHandler.h>
#include <AliNanoFilterNormalisation.h>
#include <AliVEvent.h>

#include <TChain.h>
#include <TFile.h>
#include <TObject.h>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {
  constexpr char kNames[2][16 ]{"Filter","Skimming"};
}

AliAnalysisTaskNanoAODnormalisation::AliAnalysisTaskNanoAODnormalisation(std::string taskName) : AliAnalysisTaskSE{taskName.data()},
  fOutputList{nullptr},
  fNmultBins{100},
  fMinMult{0.},
  fMaxMult{100.f},
  fCandidateEvents{nullptr},
  fSelectedEvents{nullptr},
  fCandidateEventsUE{nullptr},
  fSelectedEventsUE{nullptr}
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

  for (int iF{0}; iF < 2; ++iF) {
    fCandidateEvents[iF] = new TH2D(Form("fCandidateEvents_%s",kNames[iF]), ";Multiplicity estimator;", fNmultBins, fMinMult, fMaxMult, 5, -0.5, 4.5);
    fSelectedEvents[iF]  = (TH2D*) fCandidateEvents[iF]->Clone(Form("fSelectedEvents_%s",kNames[iF]));
    fCandidateEventsUE[iF] = new TH2D(Form("fCandidateEventsUE_%s",kNames[iF]), ";Multiplicity estimator;", fNmultBins, fMinMult, fMaxMult, 5, -0.5, 4.5);
    fSelectedEventsUE[iF]  = (TH2D*) fCandidateEvents[iF]->Clone(Form("fSelectedEventsUE_%s",kNames[iF]));

    std::string labels[5]{"Input events", "Triggered events", "Triggered + Quality cuts", "Triggered + QC + Reco vertex", "Analysis events"};
    for (int iType{0}; iType < 5; ++iType) {
      fCandidateEvents[iF]->GetYaxis()->SetBinLabel(iType + 1, labels[iType].data());
      fSelectedEvents[iF]->GetYaxis()->SetBinLabel(iType + 1, labels[iType].data());
      fCandidateEventsUE[iF]->GetYaxis()->SetBinLabel(iType + 1, labels[iType].data());
      fSelectedEventsUE[iF]->GetYaxis()->SetBinLabel(iType + 1, labels[iType].data());
    }
    fOutputList->Add(fCandidateEvents[iF]);
    fOutputList->Add(fSelectedEvents[iF]);
    fOutputList->Add(fCandidateEventsUE[iF]);
    fOutputList->Add(fSelectedEventsUE[iF]);
  }
  PostData(1, fOutputList);
}

void AliAnalysisTaskNanoAODnormalisation::FillHistograms(TH2D* candidate[2], TH2D* selected[2]) {
  AliAODInputHandler* handler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!handler) {
    ::Warning("AliAnalysisTaskNanoAODnormalisation::FillHistograms","Missing input handler");
    return;
  }
  TTree* tree = handler->GetTree();
  if (!tree) {
    ::Warning("AliAnalysisTaskNanoAODnormalisation::FillHistograms","Missing Tree");
    return;
  }

  TFile* inputFile = tree->GetCurrentFile();
  for (int iF{0}; iF < 2; ++iF) {
    AliNanoFilterNormalisation* normalisation = (AliNanoFilterNormalisation*)inputFile->Get(Form("NanoAODFilter/NanoAOD%s_scaler",kNames[iF]));
    if (!normalisation) {
      if (!iF)
        ::Fatal("AliAnalysisTaskNanoAODnormalisation::FillHistograms","Missing filtering scalers!");
      else
        ::Warning("AliAnalysisTaskNanoAODnormalisation::FillHistograms","Missing skimming scalers!");
      continue;
    }
    candidate[iF]->Add(normalisation->GetCandidateEventsHistogram());
    selected[iF]->Add(normalisation->GetSelectedEventsHistogram());
  }
  PostData(1, fOutputList); 
}

Bool_t AliAnalysisTaskNanoAODnormalisation::UserNotify() {
  FillHistograms(fCandidateEvents,fSelectedEvents);
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
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Normalisation", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(), "NanoAODNormalisation"));
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
