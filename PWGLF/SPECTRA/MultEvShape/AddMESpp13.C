#if defined(__CLING__) || defined(__ACLIC__)
#include <TTree.h>
#include <TError.h>
#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliMESbaseTask.h>
#include <AliMESpp13.h>
#include <AliMESeventInfo.h>
#endif

AliMESpp13 *AddMESpp13(Bool_t mc)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliMESpp13 *task = new AliMESpp13((char *)"MESpp13");
  mgr->AddTask(task);
  // task set-up
  task->SetMCdata(mc);

  // connect input
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // create output containers
  AliAnalysisDataContainer *co[AliMESpp13::kNcontainers] = {NULL};
  co[0] = mgr->CreateContainer("taskQA", TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  co[AliMESpp13::kTree] = mgr->CreateContainer("MES-ev", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  if (mc)
  {
    co[AliMESpp13::kMCGenTree] = mgr->CreateContainer("MES-genTrk", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
    // co[AliMESpp13::kMCMissTree] = mgr->CreateContainer("MES-missedTrk", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  }

  // connect output
  for (Int_t ios(0); ios < AliMESpp13::kNcontainers; ios++)
    if (co[ios])
      mgr->ConnectOutput(task, ios + 1, co[ios]);

  return task;
}
