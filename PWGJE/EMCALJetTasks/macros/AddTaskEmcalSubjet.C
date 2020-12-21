#if !defined(__CINT__) || defined(__MAKECINT__) || defined(__CLING__) || defined(__ROOTCLING__)
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliEmcalList.h"
#include "AliAnalysisTaskEmcalSubjet.h"
#endif

AliAnalysisTaskEmcalSubjet* AddTaskEmcalSubjet(const TString sTrks = "usedefault",
                                               const TString sClus = "usedefault",
                                               const TString sCells = "usedefault")
{

  auto mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskEmcalSubjet::AddTaskEmcalSubjet", "No analysis manager to connect to");
     return nullptr;
  }
//=============================================================================

  AliVEventHandler *pH(mgr->GetInputEventHandler());

  if (!pH) {
    ::Error("AddTaskEmcalSubjet::AddTaskEmcalSubjet", "This task requires an input event handler");
    return nullptr;
  }
//=============================================================================

  enum EDataType_t { kUnknown, kESD, kAOD };

  EDataType_t wType(kUnknown);
  if (pH->InheritsFrom("AliESDInputHandler")) wType = kESD;
  if (pH->InheritsFrom("AliAODInputHandler")) wType = kAOD;

  if (wType==kUnknown) {
    ::Error("AddTaskEmcalSubjet::AddTaskEmcalSubjet", "Unkown data input");
    return nullptr;
  }
//=============================================================================

  auto task = new AliAnalysisTaskEmcalSubjet("AliAnalysisTaskEmcalSubjet",kTRUE);
//=============================================================================

  TString sTrkName(sTrks);
  if (sTrks=="usedefault") {
    if (wType==kESD) sTrkName = "Tracks";
    if (wType==kAOD) sTrkName = "tracks";
  }

  if (sTrkName=="mcparticles") {
    task->AddMCParticleContainer(sTrkName);
  } else if ((sTrkName=="tracks") || (sTrkName=="Tracks")) {
    task->AddTrackContainer(sTrkName);
  } else if (!sTrkName.IsNull()) {
    task->AddParticleContainer(sTrkName);
  }
//=============================================================================

  TString sClsName(sClus);
  if (sClus=="usedefault") {
    if (wType==kESD) sClsName = "CaloClusters";
    if (wType==kAOD) sClsName = "caloClusters";
  }

  task->AddClusterContainer(sClsName);
//=============================================================================
  TString sCellName(sCells);
  if (sCells=="usedefault") {
    if (wType==kESD) sCellName = "EMCALCells";
    if (wType==kAOD) sCellName = "emcalCells";
  }

  task->SetCaloCellsName(sCellName);
//=============================================================================

  mgr->AddTask(task);
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("listSubjet",
                                                   AliEmcalList::Class(),
                                                   AliAnalysisManager::kOutputContainer,
                                                   AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  return task;
}
