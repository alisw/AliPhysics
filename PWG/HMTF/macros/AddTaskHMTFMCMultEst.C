#ifndef __CINT__
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskHMTFMCMultEst.h"
#endif

AliAnalysisTaskHMTFMCMultEst *AddTaskHMTFMCMultEst() {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHMTFMCMultEst", "No analysis manager to connect to.");
    return NULL;
  }

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
    "Sums",
    TList::Class(), 
    AliAnalysisManager::kOutputContainer,
    Form("%s:MultEstimators", mgr->GetCommonFileName()));

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(
    "runconditions",
    TList::Class(),
    AliAnalysisManager::kParamContainer,  // important, apparently...
    Form("%s:MultEstimators", mgr->GetCommonFileName()));
  AliAnalysisTaskHMTFMCMultEst *multEstTask = new AliAnalysisTaskHMTFMCMultEst("TaskHMTFMCMultEst");
  if (!multEstTask) {
      Error("CreateTasks", "Failed to add task!");
      return NULL;
  }
  // add estimators:
  //multEstTask->AddEstimator("Total");
  multEstTask->AddEstimator("EtaLt05");
  multEstTask->AddEstimator("EtaLt08");
  multEstTask->AddEstimator("EtaLt15");
  multEstTask->AddEstimator("Eta08_15");
  multEstTask->AddEstimator("V0A");
  multEstTask->AddEstimator("V0C");
  multEstTask->AddEstimator("V0M");

  multEstTask->SetReferenceEstimator("EtaLt05");

  // Other options
  multEstTask->SetRequireINELgt0(kFALSE);
  multEstTask->SetFillNtuple(kFALSE);

  mgr->AddTask(multEstTask);
  AliAnalysisDataContainer *inputContainer = mgr->GetCommonInputContainer();
  if(!inputContainer) {
      Error("CreateTasks", "No input container available. Failed to add task!");
      return NULL;
  }
  mgr->ConnectInput(multEstTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(multEstTask, 1, coutput1);
  mgr->ConnectOutput(multEstTask, 2, coutput2);
  return multEstTask;
}
