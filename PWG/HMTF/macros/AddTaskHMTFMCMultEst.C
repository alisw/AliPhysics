#ifndef __CINT__
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskHMTFMCMultEst.h"
#endif

AliAnalysisTaskHMTFMCMultEst *AddTaskHMTFMCMultEst(const char* globalTrigger = "") {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHMTFMCMultEst", "No analysis manager to connect to.");
    return NULL;
  }

  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("Sums%s", globalTrigger),
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 Form("%s:MultEstimators", mgr->GetCommonFileName()));

  AliAnalysisTaskHMTFMCMultEst *multEstTask = new AliAnalysisTaskHMTFMCMultEst("TaskHMTFMCMultEst");
  multEstTask->SetGlobalTrigger(globalTrigger);
  if (!multEstTask) {
      Error("CreateTasks", "Failed to add task!");
      return NULL;
  }

  mgr->AddTask(multEstTask);
  AliAnalysisDataContainer *inputContainer = mgr->GetCommonInputContainer();
  if(!inputContainer) {
      Error("CreateTasks", "No input container available. Failed to add task!");
      return NULL;
  }
  mgr->ConnectInput(multEstTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(multEstTask, 1, coutput1);

  return multEstTask;
}
