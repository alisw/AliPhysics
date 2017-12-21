#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetV0Filter.h"
#endif

AliAnalysisTaskEmcalJetV0CF *AddTaskEmcalJetV0CF(const Bool_t bH=kTRUE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskEmcalJetV0CF::AddTaskEmcalJetV0CF", "No analysis manager to connect to!!!");
     return NULL;
  }
//=============================================================================

  AliAnalysisTaskEmcalJetV0CF *task = new AliAnalysisTaskEmcalJetV0CF("AliAnaTaskEmcalJetV0CF", bH);
//=============================================================================

  mgr->AddTask(task);
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());

  if (bH) mgr->ConnectOutput(task, 1,
                             mgr->CreateContainer("listEMCalJet",
                                                  AliEmcalList::Class(),
                                                  AliAnalysisManager::kOutputContainer,
                                                  AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectOutput(task, bH ? 2 : 1,
                     mgr->CreateContainer("listEMCalJetV0CF",
                                           AliEmcalList::Class(),
                                           AliAnalysisManager::kOutputContainer,
                                           AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  return task;
}
