#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TList.h>
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSEPicoV0Maker.h"
#endif

AliAnalysisTaskSEPicoV0Maker *AddTaskPicoV0Maker(const Bool_t bMC=kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskPicoV0Maker::AddTaskPicoV0Maker", "No analysis manager to connect to!!!");
     return NULL;
  }
//=============================================================================

  AliAnalysisTaskSEPicoV0Maker *task = new AliAnalysisTaskSEPicoV0Maker("AliAnaTaskPicoV0Maker", bMC);
//=============================================================================

  mgr->AddTask(task);
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("listPicoV0MakerEH",
                                                    TList::Class(),
                                                    AliAnalysisManager::kOutputContainer,
                                                    AliAnalysisManager::GetCommonFileName()));

  if (bMC) mgr->ConnectOutput(task, 2, mgr->CreateContainer("listPicoV0MakerMC",
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  return task;
}
