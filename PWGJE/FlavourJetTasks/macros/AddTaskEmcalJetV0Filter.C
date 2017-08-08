#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include "AliEmcalList.h"
#include "AliVEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetV0Filter.h"
#endif

AliAnalysisTaskEmcalJetV0Filter *AddTaskEmcalJetV0Filter(
  const TString sTrks  = "usedefault",
  const TString sClus  = "usedefault",
  const TString sCells = "usedefault",
  const Bool_t bH = kTRUE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskEmcalJetV0Filter::AddTaskEmcalJetV0Filter", "No analysis manager to connect to!!!");
     return NULL;
  }
//=============================================================================

  AliVEventHandler *pH(mgr->GetInputEventHandler());

  if (!pH) {
    ::Error("AddTaskEmcalJetV0Filter::AddTaskEmcalJetV0Filter", "This task requires an input event handler");
    return NULL;
  }
//=============================================================================

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD };

  EDataType_t wType(kUnknown);
  if (pH->InheritsFrom("AliESDInputHandler")) wType = kESD;
  if (pH->InheritsFrom("AliAODInputHandler")) wType = kAOD;

  if (wType==kUnknown) {
    ::Error("AddTaskEmcalJetV0Filter::AddTaskEmcalJetV0Filter", "Unkown data input");
    return NULL;
  }
//=============================================================================

  TString sTrkName(sTrks);
  if (sTrks=="usedefault") {
    if (wType==kESD) sTrkName = "Tracks";
    if (wType==kAOD) sTrkName = "tracks";
  }

  TString sClsName(sClus);
  if (sClus=="usedefault") {
    if (wType==kESD) sClsName = "CaloClusters";
    if (wType==kAOD) sClsName = "caloClusters";
  }

  TString sCellName(sCells);
  if (sCells=="usedefault") {
    if (wType==kESD) sCellName = "EMCALCells";
    if (wType==kAOD) sCellName = "emcalCells";
  }
//=============================================================================

  AliAnalysisTaskEmcalJetV0Filter *task = new AliAnalysisTaskEmcalJetV0Filter("AliAnaTaskEmcalJetV0Filter", bH);

  if (sTrkName=="mcparticles") {
    task->AddMCParticleContainer(sTrkName);
  } else if ((sTrkName=="tracks") || (sTrkName=="Tracks")) {
    task->AddTrackContainer(sTrkName);
  } else if (!sTrkName.IsNull()) {
    task->AddParticleContainer(sTrkName);
  }

  task->AddClusterContainer(sClsName);
  task->SetCaloCellsName(sCellName);
//=============================================================================

  mgr->AddTask(task);
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());

  if (bH) mgr->ConnectOutput(task, 1,
                             mgr->CreateContainer("listEMCalJet",
                                                  AliEmcalList::Class(),
                                                  AliAnalysisManager::kOutputContainer,
                                                  AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  return task;
}
