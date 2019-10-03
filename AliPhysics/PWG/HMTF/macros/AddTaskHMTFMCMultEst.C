#ifndef __CINT__
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskHMTFMCMultEst.h"
#endif

AliAnalysisTaskHMTFMCMultEst *AddTaskHMTFMCMultEst(Int_t globalTrigger, const std::string &name = "") {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHMTFMCMultEst", "No analysis manager to connect to.");
    return NULL;
  }
  TString sumsName;
  if (name.size() != 0)
    sumsName = name;
  else if (globalTrigger == 0)
    sumsName = TString("SumsInel");
  else if (globalTrigger == 1)
    sumsName = TString("SumsInelGt0");
  else if (globalTrigger == 2)
    sumsName = TString("SumsV0AND");

  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(sumsName,
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
