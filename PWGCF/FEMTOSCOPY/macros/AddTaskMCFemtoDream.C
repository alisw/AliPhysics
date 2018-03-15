#include "TROOT.h"
#include "TSystem.h"

AliAnalysisTaskSE* AddTaskMCFemtoDream() {
  // the manager is static, so get the existing manager via the static method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  // just to see if all went well, check if the input event handler has been
  // connected
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }
  AliAnalysisTaskMCFemtoDream *task=
      new AliAnalysisTaskMCFemtoDream("FemtoMC");
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);
  TString QAName = Form("QA");
  AliAnalysisDataContainer *coutputQA;
  coutputQA = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      QAName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}
