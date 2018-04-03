#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliAnalysisTaskSignalLoss.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include <Rtypes.h>
#include <TString.h>
#endif

AliAnalysisTaskSignalLoss* AddTaskSignalLoss(TString tskname = "signal_loss",
TString suffix = ""){

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskSignalLoss", "No analysis manager found.");
    return 0x0;
  }
  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskSignalLoss", "This task requires an input event handler");
    return 0x0;
  }

  tskname.Append(Form("%s",suffix.Data()));

  AliAnalysisTaskSignalLoss *task = new AliAnalysisTaskSignalLoss(tskname.Data());

  float cent[14] = {-5.f,0.f,1.f,5.f,10.f,20.f,30.f,40.f,50.f,60.f,70.f,80.f,90.f,100.f};
  task->SetCentBins(13, cent);
  float pt[20] = {
    0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,1.6f,
    1.8f,2.0f,2.2f,2.6f,3.0f,3.4f,3.8f,4.4f,5.0f,6.0f
  };
  task->SetPtBins(19,pt);
  int pdgcodes[3] = {211,321,2212};
  task->SetPDGcodes(3,pdgcodes);

  mgr->AddTask(task);

  TString output = "AnalysisResults.root";
  AliAnalysisDataContainer *taskCont = mgr->CreateContainer(tskname.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      output.Data());
  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  1, taskCont);
  return task;
}
