#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TList.h>
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSEPicoV0MakerMC.h"
#endif

AliAnalysisTaskSEPicoV0MakerMC *AddTaskPicoV0MakerMC()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskPicoV0MakerMC::AddTaskPicoV0MakerMC", "No analysis manager to connect to!!!");
     return NULL;
  }
//=============================================================================

  AliAnalysisTaskSEPicoV0MakerMC *task = new AliAnalysisTaskSEPicoV0MakerMC("AliAnaTaskPicoV0MakerMC");
//=============================================================================

  mgr->AddTask(task);
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("listPicoV0MakerMCdet",
                                                    TList::Class(),
                                                    AliAnalysisManager::kOutputContainer,
                                                    AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  return task;
}
