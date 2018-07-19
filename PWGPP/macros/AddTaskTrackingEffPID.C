#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskTrackingEffPID.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliPID.h"
#endif

AliAnalysisTaskTrackingEffPID* AddTaskTrackingEffPID(TString suffix = "",
						     TString collSyst="pp",
						     bool useGeneratedKine=kTRUE) {

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskLFefficiencies", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskLFefficiencies", "This task requires an input event handler");
    return 0x0;
  }

  AliAnalysisTaskTrackingEffPID *taskeff = new AliAnalysisTaskTrackingEffPID();
  taskeff->SetUseGeneratedKine(useGeneratedKine);
  taskeff->SetCollisionSystem(collSyst);

  mgr->AddTask(taskeff);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":TrackEffPID";

  TString listname="listTrackEffPID";
  listname+=suffix.Data();

  AliAnalysisDataContainer *coutput = mgr->CreateContainer(listname,
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   outputFileName.Data() );

  mgr->ConnectInput  (taskeff,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (taskeff,  1, coutput);
  return taskeff;
}
