#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSEPicoV0Filter.h"
#endif

AliAnalysisTaskSEPicoV0Filter *AddTaskPicoV0Filter(const Bool_t bMC=kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskPicoV0Filter::AddTaskPicoV0Filter", "No analysis manager to connect to!!!");
     return NULL;
  }
//=============================================================================

  AliAnalysisTaskSEPicoV0Filter *task = new AliAnalysisTaskSEPicoV0Filter("AliAnaTaskPicoV0Filter");
  if (bMC) task->SetAnaInfoMC(bMC);
//=============================================================================

  mgr->AddTask(task);
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
/*mgr->ConnectOutput(task, 1, mgr->CreateContainer("listPicoV0Filter",
                                                    TList::Class(),
                                                    AliAnalysisManager::kOutputContainer,
                                                    AliAnalysisManager::GetCommonFileName()));*/
//=============================================================================

  return task;
}
