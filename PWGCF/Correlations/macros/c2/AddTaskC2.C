#ifndef __CINT__
#include "TDirectory.h"
#include "TList.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskC2.h"
#include "AliAnalysisTaskC2Base.h"
#include "AliVEvent.h"

#include <iostream>
#include <string>
#endif

AliAnalysisTaskC2 *AddTaskC2(Int_t mode) {
  TString postfix;
  if (mode == AliAnalysisTaskC2Base::kMCTRUTH){
    ::Info("AddTaskC2", "Running in MC Truth mode");
    postfix = "_mcTruth";
  }
  else if (mode == AliAnalysisTaskC2::kRECON){
    ::Info("AddTaskC2", "Running in Reconstructed mode");
    postfix = "_recon";
  }
  else {
    ::Error("AddTaskC2", "Invalid mode");
    return NULL;
  }

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskC2", "No analysis manager to connect to.");
    return NULL;
  }
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("sums%s", postfix.Data()) ,
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 Form("%s:c2_correlations%s",
			      mgr->GetCommonFileName(),
			      postfix.Data()));

  AliAnalysisTaskC2 *c2Task = new AliAnalysisTaskC2("TaskC2", mode);
  //c2Task->SelectCollisionCandidates(AliVEvent::kAny);

  if (!c2Task) {
      Error("CreateTasks", "Failed to add task!");
      return NULL;
  }

  mgr->AddTask(c2Task);
  AliAnalysisDataContainer *inputContainer = mgr->GetCommonInputContainer();
  if(!inputContainer) {
      Error("CreateTasks", "No input container available. Failed to add task!");
      return NULL;
  }
  mgr->ConnectInput(c2Task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(c2Task, 1, coutput1);

  return c2Task;
}
