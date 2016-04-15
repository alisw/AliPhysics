#ifndef __CINT__
#include "TDirectory.h"
#include "TList.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskC2.h"
#include "AliAnalysisTaskC2Base.h"
#include "AliAnalysisC2Settings.h"
#include "AliVEvent.h"

#include <iostream>
#include <string>
#endif

AliAnalysisTaskC2 *AddTaskC2() {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskC2", "No analysis manager to connect to.");
    return NULL;
  }
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("sums") ,
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 Form("%s:c2_correlations", mgr->GetCommonFileName()));

  AliAnalysisTaskC2 *c2Task = new AliAnalysisTaskC2("TaskC2");
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
