#include "AliAnalysisTaskNanoAODskimming.h"
#include <AliAnalysisManager.h>
#include <AliAODExtension.h>
#include <AliAODHandler.h>
#include <AliESDtrackCuts.h>
#include <AliInputEventHandler.h>
#include <AliNanoFilterNormalisation.h>
#include <AliVEvent.h>

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
  DefineOutput(1, AliNanoFilterNormalisation::Class());
}

AliAnalysisTaskNanoAODskimming::~AliAnalysisTaskNanoAODskimming() {
  if (fNormalisation)
    delete fNormalisation;
}

void AliAnalysisTaskNanoAODskimming::UserCreateOutputObjects() {
  std::string normName = std::string(fName) + "_scaler";
  fNormalisation = new AliNanoFilterNormalisation(normName.data(), normName.data());
  PostData(1, fNormalisation);
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

  int anyBin = fNormalisation->GetCandidateEventsHistogram()->GetYaxis()->FindBin(AliNanoFilterNormalisation::kAnyEvent);
  double nCan = fNormalisation->GetCandidateEventsHistogram()->Integral(0,-1,anyBin,anyBin);
  int analysisBin = fNormalisation->GetSelectedEventsHistogram()->GetYaxis()->FindBin(AliNanoFilterNormalisation::kAnalysisEvent);
  double nSel = fNormalisation->GetSelectedEventsHistogram()->Integral(0,-1,analysisBin,analysisBin);
  std::cout << "****************************************************************" << std::endl;
  std::cout << "* AliAnalysisTaskNanoAODskimming summary" << std::endl;
  std::cout << "* Number of processed events:" << nCan << std::endl;
  std::cout << "* Number of selected events:" << nSel << std::endl;
  std::cout << "****************************************************************" << std::endl;

  AliAODHandler* handler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  if (!handler) return;
  AliAODExtension *extNanoAOD = handler->GetFilteredAOD("AliAOD.NanoAOD.root");
  if (!extNanoAOD) return;
  TTree* nanoTree = extNanoAOD->GetTree();
  if (!nanoTree) return;
  TList* userInfo = nanoTree->GetUserInfo();
  if (!userInfo) return;
  userInfo->Add(fNormalisation->Clone());
}

void AliAnalysisTaskNanoAODskimming::Terminate(Option_t *) {

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
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Skimming_Normalisation", AliNanoFilterNormalisation::Class(), AliAnalysisManager::kOutputContainer, "AliAOD.NanoAOD.root:NanoAODskimming");
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
