#include "AliAnalysisTaskNanoAODskimming.h"
#include <AliAnalysisManager.h>
#include <AliAODExtension.h>
#include <AliAODHandler.h>
#include <AliESDtrackCuts.h>
#include <AliInputEventHandler.h>
#include <AliNanoFilterNormalisation.h>
#include <AliVEvent.h>
#include <AliVTrack.h>

#include <TChain.h>
#include <TObject.h>

#include <algorithm>
#include <cmath>
#include <iostream>

AliAnalysisTaskNanoAODskimming::AliAnalysisTaskNanoAODskimming(std::string taskName) : AliAnalysisTaskSE{taskName.data()},
  fEventCuts{},
  fNormalisation{nullptr},
  fOtherEventCuts{}
{
  fOtherEventCuts.SetOwner(true);
  DefineInput(0, TChain::Class());  
}

AliAnalysisTaskNanoAODskimming::~AliAnalysisTaskNanoAODskimming() {
  delete fNormalisation;
}

void AliAnalysisTaskNanoAODskimming::UserCreateOutputObjects() {
  std::string normName = std::string(fName) + "_scaler";
  fNormalisation = new AliNanoFilterNormalisation(normName.data(), normName.data());
}

void AliAnalysisTaskNanoAODskimming::UserExec(Option_t* /*option*/) {
  AliVEvent *ev = InputEvent();

  bool acceptEvent = fEventCuts.IsSelected(ev);
  bool selMask[4]{
    fEventCuts.GetAliEventCuts().CheckNormalisationMask(AliEventCuts::kTriggeredEvent),
    fEventCuts.GetAliEventCuts().CheckNormalisationMask(AliEventCuts::kPassesNonVertexRelatedSelections),
    fEventCuts.GetAliEventCuts().CheckNormalisationMask(AliEventCuts::kHasReconstructedVertex),
    fEventCuts.GetAliEventCuts().CheckNormalisationMask(AliEventCuts::kPassesAllCuts)
  };
  float v0m = fEventCuts.GetAliEventCuts().GetCentrality() < 0 ? 0.5f : fEventCuts.GetAliEventCuts().GetCentrality();
  fNormalisation->FillCandidate(selMask[0], selMask[1], selMask[2], selMask[3], v0m);

  TIter cutsIter(&fOtherEventCuts);
  AliAnalysisCuts* cut = nullptr;
  while ((cut = (AliAnalysisCuts*)cutsIter()) && acceptEvent) {
    acceptEvent = cut->IsSelected(ev);
  }

  if (!acceptEvent)
    AliAnalysisManager::GetAnalysisManager()->BreakExecutionChain();
  else
    fNormalisation->FillSelected(selMask[0], selMask[1], selMask[2], selMask[3], v0m);
  
}

void AliAnalysisTaskNanoAODskimming::FinishTaskOutput() {
  // We save here the user info

  AliAODHandler* handler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  AliAODExtension *extNanoAOD = handler->GetFilteredAOD("AliAOD.NanoAOD.root");
  extNanoAOD->GetTree()->GetUserInfo()->Add(fNormalisation);
}

AliAnalysisTaskNanoAODskimming* AliAnalysisTaskNanoAODskimming::AddTask(std::string name) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return nullptr;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPhysicsSelection", "This task requires an input event handler");
    return nullptr;
  }

  AliAnalysisTaskNanoAODskimming *task = new AliAnalysisTaskNanoAODskimming(name);
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  return task;
}